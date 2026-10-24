r"""Gates of the diffusion realizer table ``DIFFUSION_REALIZERS`` and of ``ZeroFluxCondition`` (wave 2, WF2.3 slice S3).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.7 ("Diffusion: ``zero_flux.py`` and ``DIFFUSION_REALIZERS``"), section
4.2 (the diffusion column), section 3 (law L-L: the realization is a functor on the response network) and section 6.3 (the
eight done-when gates; this file holds gates 2, 3, 4, 6 and 7, and the structural rows the table owes). Gate 1 (the pin rows) is
in ``tests/gates/transport/test_admission_matrix.py``, gate 5 (the carry) in ``test_diffusion_boundary_carry.py`` and gate 8 (the
qualname scan) in ``test_condition_family_qualnames.py``.

The claims, all read through the door a problem uses, ``DIFFUSION_REALIZERS.admit`` and ``.realize(declared, site, space)``, on
the metadata-free method space (the albedo-family realizers read nothing from it), plus the mesh-level legs that read the same
table through ``DiffusionMesh``:

* **Diffusion times the empty network is TIER 1** (``Realized(faithful=True)``): the deciding line of the retired realizer was
  ``return response.amplitude`` (``measure_diffusion_vacuum.md``): the vacuum law has amplitude 0, ``J^- = 0`` maps to ``A = 0``
  exactly, and a scalar commutes with every linear projection. The Marshak closure is the METHOD's (the P1 closure applied
  identically to every albedo), not the vacuum law's admission. The specular and the diffuse chains are TIER 2
  (``Realized(faithful=False)``): diffusion annihilates the shape, so ``specular(0.4)`` and ``white(0.4)`` have one image. Zero flux
  is tier 1 (``A = -1``) and diffusion-only.
* **A mixed wall is the functor's sum of its chains' scalars** (``A = s + d``): the image of a sum is the sum of the images, the
  chains summed in ``CHAIN_ORDER``. The amplitudes of the row have ``s + d != 1``, because the mutation "rescale the diffuse image
  by the total return fraction" is inert at ``rho = 1``.
* **A source is a plumbing refusal** (``Refused(PLUMBING, #656)``): the fixed-source arm lifts with no theory change. It is
  refused BEFORE sampling, at every site, alone and over a network.
* **The translation deck is a topology refusal**: a per-face scalar has no slot for cross-face coupling; the mirror deck is
  tier 2 (``A = 1``).
* **The albedo above one is unspellable** (#265's albedo item): the cone half is the S2a package's (green today), so this
  gate's own red input is the other half: the old door that realized ``albedo = 5.0`` (``BOUNDARY_OPERATOR_REGISTRY``,
  ``DiffusionBoundaryRealizer``, ``_partial_current_albedo``) is gone, and no admitted wall's ``k`` exceeds the mirror's.
* **The unregistered door** of a signed-current table cites #659 and not the state cone: a negative return fraction is
  inside diffusion's cone.

First reds, at the base tip (``0757a8d9``): ``ImportError`` from ``import orpheus.diffusion.*`` (the package ``__init__`` imports
``augmented_mesh``, which imports ``resolve_boundary_conditions`` from ``orpheus.transport.method``, removed by S2b), then
``ModuleNotFoundError: orpheus.diffusion.zero_flux`` and the absent name ``DIFFUSION_REALIZERS``. The one row that is green at
the base is the cone half of ``test_super_markov_albedo_is_unspellable``'s first leg, which is S2a's and is the positive control of
the leg that is red.

The ladder::

    test_the_table_registers_the_six_realizers_and_is_sealed      foundation: the six realizer classes, the sealed table
    test_the_unregistered_door_cites_the_tables_issue_on_the_signed_cone   foundation: #659, no STATE_CONE
        ^
    test_diffusion_vacuum_is_tier_1                               edge: the empty network and zero flux, tier 1 and tier 2 chains
    test_diffusion_translation_is_topological                     edge: the decks
    test_diffusion_source_is_plumbing                             edge: the source refusal
        ^
    test_a_mixed_wall_realizes_as_the_sum_of_scalars              interior: law L-L on the diffusion table
    test_super_markov_albedo_is_unspellable                       interior: through the solve
    test_the_mesh_reads_the_binding_through_the_shared_door       composition: DiffusionMesh, the shared realize_boundary

Mutations for the reviewers (each on a COPY of ``boundary_realizer.py`` / ``zero_flux.py``, the verdict per arm):

* M-vacuum-tier-2: ``DiffusionZero.status`` and its admit return ``Realized(False)`` (the Marshak ``Approximated`` of DRAFT 1): reddens
  ``test_diffusion_vacuum_is_tier_1`` and the vacuum row of the pin.
* M-chain-tier-1: ``DiffusionSpecularChain.admit`` returns ``Realized(True)`` (the control named in the measurement): reddens
  ``test_diffusion_vacuum_is_tier_1`` (specular row) and the specular rows of the pin.
* M-rescale (law L-L): the diffuse chain's image is rescaled by ``s + d`` (read from the weight of the other chain): inert at
  ``s + d = 1``, reddens ``test_a_mixed_wall_realizes_as_the_sum_of_scalars`` at ``(0.3, 0.5)``.
* M-first-chain-only: the functor realizes only the first chain of a mixed wall (or sums in the order the wall was built):
  reddens the sum row and the order row.
* M-source-admitted: ``DiffusionSource.admit`` returns ``Realized(True)``: reddens ``test_diffusion_source_is_plumbing`` (every leg)
  and the inflow rows of the pin. M-source-cites-290: the issue 656 replaced by the closed 290: reddens the same.
* M-translation-admitted: ``DiffusionDecks.admit`` returns ``Realized(False)`` for ``TranslationDeck``: reddens
  ``test_diffusion_translation_is_topological`` (both the table leg and the mesh leg) and the periodic-slab construction.
* M-zero-flux-positive: ``DiffusionZeroFlux`` realizes ``+IdentityOperator()``: reddens the zero-flux action leg and the carry rows
  ``keff_diffusion/*/zero_flux``.
* M-door-cone: the table built with ``state_cone=StateCone.NONNEGATIVE_FLUX``: reddens ``test_the_table_registers_...`` and
  ``test_the_unregistered_door_...`` (a negative probe becomes ``STATE_CONE``).
* M-old-door-survives: ``BOUNDARY_OPERATOR_REGISTRY`` (or ``_partial_current_albedo``) kept beside the table: reddens
  ``test_super_markov_albedo_is_unspellable``'s retirement leg.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tests.gates.diffusion import _wave2_helpers as H
from tests.gates.transport._boundary_conditions_helpers import local

_THIS = "tests/gates/diffusion/test_diffusion_realizer_table.py"
_TABLE = f"{_THIS}::test_the_table_registers_the_six_realizers_and_is_sealed"


def _realize(declared: Any, kind: str = "CART") -> Any:
    return H.table().realize(declared, H.site(kind), H.minimal_space())


# ═════════════════════════════════════════════════════════════════════════════
# The table itself
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_the_table_registers_the_six_realizers_and_is_sealed() -> None:
    r"""``DIFFUSION_REALIZERS`` is a sealed ``RealizerTable`` of method ``"diffusion"`` on the signed-current cone whose
    unregistered-unit door cites #659; it is filled by ``DiffusionZero`` (``ZeroNetwork``), ``DiffusionSpecularChain``,
    ``DiffusionDiffuseChain``, ``DiffusionZeroFlux`` (``ZeroFluxCondition``), ``DiffusionSource`` (a ``SourceRealizer`` that
    overrides no ``sample``: a source is an array and the door raises before the default body) and ``DiffusionDecks`` (the one
    ``DeckRealizer``). The catalogue lists the zero network and zero flux at ``Realized(True)``, the two chains at
    ``Realized(False)`` and the source at ``Refused(PLUMBING, 656)``. The retired door is gone from the realizer module:
    ``DiffusionBoundaryRealizer`` and ``_partial_current_albedo``; ``_albedo_operator`` stays.

    Red input: the base tip has no table (``import orpheus.diffusion`` raises, and the module defines
    ``DiffusionBoundaryRealizer``). Mutation: M-door-cone; a table sealed before ``DiffusionDecks`` registers; ``DiffusionSource``
    given a ``sample``.
    """
    module = H.diffusion_module("boundary_realizer")
    B = H.bc()
    table = H.table()
    assert table.method == "diffusion" and table.state_cone is B.StateCone.SIGNED_CURRENT and table.unregistered_issue == 659
    for name, base in (
        ("DiffusionZero", B.ResponseRealizer), ("DiffusionSpecularChain", B.ResponseRealizer), ("DiffusionDiffuseChain", B.ResponseRealizer),
        ("DiffusionZeroFlux", B.ResponseRealizer), ("DiffusionSource", B.SourceRealizer), ("DiffusionDecks", B.DeckRealizer),
    ):
        assert issubclass(getattr(module, name), base), f"{name} is not a {base.__name__}"
    assert "sample" not in vars(module.DiffusionSource) and not hasattr(module.DiffusionSource, "realize"), "a source realizer has no body to stub"
    assert not hasattr(module, "DiffusionBoundaryRealizer") and not hasattr(module, "_partial_current_albedo"), "the retired door must be gone"
    assert callable(module._albedo_operator), "_albedo_operator stays (the shared scalar-to-operator collapse)"
    catalogue = table.catalogue()
    assert H.summary(catalogue[B.ZeroNetwork]) == ("R", True) and H.summary(catalogue[H.zero_flux_condition()]) == ("R", True)
    assert H.summary(catalogue[B.SpecularChain]) == ("R", False) and H.summary(catalogue[B.DiffuseChain]) == ("R", False)
    assert H.summary(catalogue[B.InflowTerms]) == ("F", "PLUMBING", 656)


@pytest.mark.foundation
@pytest.mark.rests_on(_TABLE)
def test_the_unregistered_door_cites_the_tables_issue_on_the_signed_cone() -> None:
    r"""A condition with no entry in the diffusion table is ``Refused(PLUMBING, #659)`` whatever its return fraction: a negative
    return fraction (SN's ``STATE_CONE`` case) is INSIDE the signed-current cone of diffusion, so the state-cone axis never fires
    here. The test-local ``NegativeProbe`` (``return_fraction == -1``, no entry) stands for a method-only condition of another
    method, and ``Probe(0.5)`` for a non-negative one; ``require`` raises ``BoundaryRefused`` carrying the plumbing admission.

    Red input: the base tip has no table. Mutation M-door-cone (the table on the non-negative cone reads ``STATE_CONE``).
    """
    B = H.bc()
    probes = local()
    where = H.site("CART")
    for probe in (probes.NegativeProbe(), probes.Probe(0.5)):
        admission = H.table().admit(probe, where)
        assert H.summary(admission) == ("F", "PLUMBING", 659), f"{type(probe).__name__}: {admission!r}"
    with pytest.raises(B.BoundaryRefused) as info:
        H.table().require(probes.NegativeProbe(), where)
    assert info.value.admission.axis is B.RefusalAxis.PLUMBING and info.value.admission.issue == 659


# ═════════════════════════════════════════════════════════════════════════════
# Gate 2: the tiers
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_diffusion_vacuum_is_tier_1() -> None:
    r"""Gate 2 of 6.3. ``WallCondition.vacuum()`` (and the two other spellings of the empty network, ``specular(0)`` and
    ``white(0)``) realizes ``ZeroOperator`` with admission ``Realized(True)`` at every site with a trace (``H.TRACE_SITES``; the
    hollow inner site's refusals are the admission matrix's diffusion column); ``specular(0.4)`` and ``white(0.4)``
    realize the same ``0.4 * IdentityOperator()`` (the same action on the probe, bit for bit, equal to ``0.4 * v``) with admission
    ``Realized(False)``; ``ZeroFluxCondition()`` realizes ``-IdentityOperator()`` (the action ``-v``) with ``Realized(True)``.
    Every realized operator carries ``BlockRole.BOUNDARY``. The pairs ``(admit, realize)`` are read from the same call, so a
    realizer whose body and admission disagree is seen.

    Red input: the base tip states none of these admissions (the table does not exist; the retired realizer returned a bare
    operator with no admission). Mutations: M-vacuum-tier-2, M-chain-tier-1, M-zero-flux-positive.
    """
    from orpheus.numerics.operator import BlockRole, ZeroOperator

    B = H.bc()
    W = B.WallCondition
    zero_flux = H.zero_flux_condition()()
    for kind in H.TRACE_SITES:
        where = H.site(kind)
        for label, declared in (("vacuum", W.vacuum()), ("specular(0)", W.specular(0.0)), ("white(0)", W.white(0.0))):
            assert H.summary(H.table().admit(declared, where)) == ("R", True), f"{label} @ {kind}"
            op = H.table().realize(declared, where, H.minimal_space())
            assert isinstance(op, ZeroOperator) and op.block_role is BlockRole.BOUNDARY, f"{label} @ {kind}: {type(op).__name__}"
            np.testing.assert_array_equal(H.realized_action(op), 0.0)
        images = []
        for label, declared in (("specular(0.4)", W.specular(0.4)), ("white(0.4)", W.white(0.4))):
            assert H.summary(H.table().admit(declared, where)) == ("R", False), f"{label} @ {kind}"
            op = H.table().realize(declared, where, H.minimal_space())
            assert op.block_role is BlockRole.BOUNDARY
            np.testing.assert_array_equal(H.realized_action(op), 0.4 * H.PROBE)
            images.append(H.realized_action(op))
        np.testing.assert_array_equal(images[0], images[1])  # the tier-2 witness: one image for two spellings
        assert H.summary(H.table().admit(zero_flux, where)) == ("R", True), f"zero flux @ {kind}"
        op = H.table().realize(zero_flux, where, H.minimal_space())
        assert op.block_role is BlockRole.BOUNDARY
        np.testing.assert_array_equal(H.realized_action(op), -H.PROBE)


# ═════════════════════════════════════════════════════════════════════════════
# Gate 6: law L-L
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_a_mixed_wall_realizes_as_the_sum_of_scalars() -> None:
    r"""Gate 6 of 6.3 (law L-L on the diffusion table). The image of a mixed wall is the sum of the images of its chains:
    ``from_amplitudes(s, d)`` applies as ``s v + d v``, i.e. ``A = s + d``, admission ``Realized(False)`` (the meet of two tier-2
    units). Each single-chain image is ``array_equal`` to ``w v`` (a scalar times a vector reorders no addition); the sum is
    ``nulp``-close to the sum of the two single images and to ``(s + d) v`` (one added pair: depth 1, so 2 ulp), and is built
    in ``CHAIN_ORDER``, so ``specular(s) + white(d)`` and ``white(d) + specular(s)`` realize the SAME action bit for bit. The
    rows ``(0.3, 0.5)`` and ``(0.7, 0.1)`` have ``s + d != 1``, with ``s != d``, so a rescale by the total return fraction (M-rescale)
    and a swap of the two weights are both seen; ``(0.3, 0.7)`` (``s + d = 1``, the identity as a sum) is the control that is
    blind to a rescale.

    Red input: the base tip has no table. Mutations: M-rescale (inert on the third row, red on the first two), M-first-chain-only.
    """
    B = H.bc()
    W = B.WallCondition
    where = H.site("CART")
    for s, d in ((0.3, 0.5), (0.7, 0.1), (0.3, 0.7)):
        mixed = W.from_amplitudes(s, d)
        assert H.summary(H.table().admit(mixed, where)) == ("R", False), f"({s}, {d})"
        image = H.realized_action(H.table().realize(mixed, where, H.minimal_space()))
        specular = H.realized_action(H.table().realize(W.specular(s), where, H.minimal_space()))
        diffuse = H.realized_action(H.table().realize(W.white(d), where, H.minimal_space()))
        np.testing.assert_array_equal(specular, s * H.PROBE)
        np.testing.assert_array_equal(diffuse, d * H.PROBE)
        np.testing.assert_array_almost_equal_nulp(image, specular + diffuse, nulp=2)
        np.testing.assert_array_almost_equal_nulp(image, (s + d) * H.PROBE, nulp=2)
        built = {"specular first": W.specular(s) + W.white(d), "diffuse first": W.white(d) + W.specular(s)}
        for label, wall in built.items():
            np.testing.assert_array_equal(H.realized_action(H.table().realize(wall, where, H.minimal_space())), image, err_msg=label)


# ═════════════════════════════════════════════════════════════════════════════
# Gate 7: the source
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_diffusion_source_is_plumbing(monkeypatch: pytest.MonkeyPatch) -> None:
    r"""Gate 7 of 6.3. A nonzero source is ``Refused(PLUMBING, 656)`` at every site with a trace (the hollow inner site
    refuses it on #511 first, the admission matrix's row), alone (``WallCondition.uniform(1)``, the bare
    ``InflowTerms``), over a response (``specular(0.5) + uniform(1)``) and through ``realize``; ``sample_inflow`` raises
    ``BoundaryRefused`` BEFORE any sampling (``sample_terms`` is called 0 times; a positive control calls the default
    ``SourceRealizer.sample`` of the registered realizer directly and counts 1, so the spy is live), and a ``DiffusionMesh`` over a
    uniform inflow refuses at construction with the same admission. #656 is the issue that tracks diffusion's affine channel
    (the closed #290 is not it).

    Red input: the base tip has no table and the retired realizer raised ``BoundaryError`` for ``PrescribedInflow`` with no issue.
    Mutations: M-source-admitted, M-source-cites-290.
    """
    B = H.bc()
    W = B.WallCondition
    plumbing = ("F", "PLUMBING", 656)
    for kind in H.TRACE_SITES:
        where = H.site(kind)
        for label, declared in (("uniform", W.uniform(1.0)), ("terms", B.InflowTerms(1.0, ())), ("response over source", W.specular(0.5) + W.uniform(1.0))):
            assert H.summary(H.table().admit(declared, where)) == plumbing, f"{label} @ {kind}"
        with pytest.raises(B.BoundaryRefused) as info:
            H.table().realize(W.uniform(1.0), where, H.minimal_space())
        assert H.summary(info.value.admission) == plumbing

    import orpheus.transport.boundary_conditions.realization as realization

    calls: list[Any] = []
    real = realization.sample_terms
    monkeypatch.setattr(realization, "sample_terms", lambda terms, space: calls.append(terms) or real(terms, space))
    terms = B.InflowTerms(1.0, ())
    space = type("Space", (), {"shape": (2,)})()
    with pytest.raises(B.BoundaryRefused):
        H.table().sample_inflow(terms, H.site("CART"), space)
    assert calls == [], "the source verb sampled before the door refused"
    registered = H.diffusion_module("boundary_realizer").DiffusionSource()
    np.testing.assert_array_equal(B.SourceRealizer.sample(registered, terms, H.site("CART"), space), 1.0)
    assert len(calls) == 1, "the spy is dead: the default body did not reach sample_terms"

    mesh = H.chart_mesh("slab", left=0, right=0)
    with pytest.raises(B.BoundaryRefused) as info:
        H.diffusion_mesh(mesh, {0: W.uniform(1.0)})
    assert H.summary(info.value.admission) == plumbing


# ═════════════════════════════════════════════════════════════════════════════
# Gate 4: the decks
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_diffusion_translation_is_topological() -> None:
    r"""Gate 4 of 6.3. ``TranslationDeck`` is ``Refused(TOPOLOGY)`` (a per-face scalar has no slot for cross-face coupling) at
    a Cartesian site, ``realize`` raises ``BoundaryRefused`` carrying that admission, and a slab with ``TranslationDeck`` on both
    faces refuses at ``DiffusionMesh`` construction (the periodic slab, which the retired realizer refused with a bare
    ``BoundaryError``). ``MirrorDeck`` is ``Realized(False)`` (tier 2) and realizes ``IdentityOperator()`` (``A = 1``: the
    action ``v``, the type ``IdentityOperator``), stamped ``BlockRole.BOUNDARY``.

    Red input: the base tip has no table. Mutation: M-translation-admitted.
    """
    from orpheus.geometry.deck import MirrorDeck, TranslationDeck
    from orpheus.numerics.operator import BlockRole, IdentityOperator

    B = H.bc()
    where = H.site("CART")
    refusal = H.table().admit(TranslationDeck(), where)
    assert H.summary(refusal) == ("F", "TOPOLOGY", None), f"{refusal!r}"
    with pytest.raises(B.BoundaryRefused) as info:
        H.table().realize(TranslationDeck(), where, H.minimal_space())
    assert H.summary(info.value.admission) == ("F", "TOPOLOGY", None)
    assert H.summary(H.table().admit(MirrorDeck(), where)) == ("R", False)
    op = H.table().realize(MirrorDeck(), where, H.minimal_space())
    assert type(op) is IdentityOperator and op.block_role is BlockRole.BOUNDARY
    np.testing.assert_array_equal(H.realized_action(op), H.PROBE)

    with pytest.raises(B.BoundaryRefused) as info:
        H.diffusion_mesh(H.chart_mesh("slab", left=TranslationDeck(), right=TranslationDeck()))
    assert info.value.admission.axis is B.RefusalAxis.TOPOLOGY


# ═════════════════════════════════════════════════════════════════════════════
# Gate 3: the albedo above one
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_super_markov_albedo_is_unspellable() -> None:
    r"""Gate 3 of 6.3 (#265's albedo item). Two legs, of which only the second is red at the base.

    *The cone leg* (S2a's, green today): ``WallCondition.specular(5.0)``, ``.white(5.0)`` and ``from_amplitudes(3.0, 2.0)`` raise
    ``SubmarkovViolationError`` at construction, and ``white(0.6) + white(0.6)`` is a ``WallDisplacement`` that
    ``DiffusionMesh`` refuses to bind (``TypeError``), so no spelling of an albedo above 1 reaches the realizer. It is this
    leg's positive control that the same constructions with ``0.6`` (and ``1.0``) build.

    *The retirement-and-bound leg*: the door through which the old code realized ``albedo = 5.0`` (``k = 2.972`` on the slab at
    ``5d04c296``) is gone (no ``DiffusionMesh.BOUNDARY_OPERATOR_REGISTRY``, no ``DiffusionBoundaryRealizer``, no
    ``_partial_current_albedo``), and through the solve every admitted wall's ``k`` is bounded by the mirror's: for
    ``white(a)`` with ``a`` in ``(0.3, 0.7, 1.0)`` the ``k`` is strictly increasing and ``k(white(1.0))`` equals the mirror's
    (a tier-2 image: the same operator), on the two-region two-group slab.

    Red input: the base tip raises ``ImportError`` at the first diffusion read. Mutation: M-old-door-survives; an albedo
    realized unclamped (a ``WallCondition`` that admits ``specular(5.0)``) reddens the cone leg.
    """
    from orpheus.geometry.deck import MirrorDeck

    B = H.bc()
    W = B.WallCondition
    for build in (lambda: W.specular(5.0), lambda: W.white(5.0), lambda: W.from_amplitudes(3.0, 2.0)):
        with pytest.raises(B.SubmarkovViolationError):
            build()
    for control in (W.specular(0.6), W.white(0.6), W.white(1.0), W.from_amplitudes(0.6, 0.4)):
        assert isinstance(control, W)
    displacement = W.white(0.6) + W.white(0.6)
    assert isinstance(displacement, B.WallDisplacement)
    with pytest.raises(TypeError):
        H.diffusion_mesh(H.chart_mesh("slab", left=0, right=0), {0: displacement})  # type: ignore[dict-item]

    augmented = H.diffusion_module("augmented_mesh")
    assert not hasattr(augmented.DiffusionMesh, "BOUNDARY_OPERATOR_REGISTRY"), "the string-tag registry is the old door of albedo = 5.0"
    assert not hasattr(H.diffusion_module("boundary_realizer"), "DiffusionBoundaryRealizer")

    slab = H.chart_mesh("slab", left=0, right=0)
    ks = [H.solve(slab, {0: W.white(a)}).keff for a in (0.3, 0.7, 1.0)]
    mirror = H.solve(H.chart_mesh("slab", left=MirrorDeck(), right=MirrorDeck())).keff
    assert ks == sorted(ks) and len(set(ks)) == 3, f"k(white(a)) is not strictly increasing in a: {ks}"
    assert ks[-1] == mirror, f"white(1.0) and the mirror deck are one image in diffusion: {ks[-1]!r} != {mirror!r}"


# ═════════════════════════════════════════════════════════════════════════════
# The mesh reads the table through the shared door
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.rests_on(_TABLE)
def test_the_mesh_reads_the_binding_through_the_shared_door(monkeypatch: pytest.MonkeyPatch) -> None:
    r"""``DiffusionMesh`` is a ``TransportMethod`` as ``SNProblem`` is: it holds ``boundary_realizers is DIFFUSION_REALIZERS``, the
    ``binding`` of its ``boundary_conditions`` (``bind(carrier, conditions)``) and a ``bc`` realized through the shared
    ``orpheus.transport.method.realize_boundary`` (a counting spy on the name the module imports: exactly one call per
    construction, PATH not output), whose ``realize_boundary_role`` arm realizes on the face's method space. A deck face and a
    tagged face are both in ``bc``; the tagged face's operator acts as the declared scalar; the old resolver is gone from
    ``orpheus.transport.method`` and from the mesh module. A face tagged with no condition declared is ``UnboundTagError`` (no
    default wall).

    Red input: the base tip's mesh module imports ``resolve_boundary_conditions``. Mutation: ``bc`` built by a diffusion-local
    loop (the spy reads 0 calls); ``boundary_realizers`` set to a table copy.
    """
    from orpheus.geometry.deck import MirrorDeck
    from orpheus.transport.method import TransportMethod

    B = H.bc()
    W = B.WallCondition
    augmented = H.diffusion_module("augmented_mesh")
    calls: list[Any] = []
    real = augmented.realize_boundary
    monkeypatch.setattr(augmented, "realize_boundary", lambda method: calls.append(method) or real(method))
    mesh = H.diffusion_mesh(H.chart_mesh("slab", left=MirrorDeck(), right=0), {0: W.white(0.5)})
    assert len(calls) == 1 and calls[0] is mesh, "the mesh did not realize its binding through the shared realize_boundary exactly once"
    assert isinstance(mesh, TransportMethod)
    assert mesh.boundary_realizers is H.table() and isinstance(mesh.binding, B.Binding)
    assert set(mesh.bc) == {"xmin", "xmax"} and set(mesh.binding.decks) == {"xmin"} and dict(mesh.binding.tags) == {"xmax": 0}
    np.testing.assert_array_equal(H.realized_action(mesh.bc["xmin"]), H.PROBE)  # the mirror deck, A = 1
    np.testing.assert_array_equal(H.realized_action(mesh.bc["xmax"]), 0.5 * H.PROBE)  # the white(0.5) tag, A = 0.5
    arm = mesh.realize_boundary_role(W.white(0.5), "xmax")
    np.testing.assert_array_equal(H.realized_action(arm), 0.5 * H.PROBE)

    method = H.diffusion_module("augmented_mesh")
    import orpheus.transport.method as shared

    assert not hasattr(shared, "resolve_boundary_conditions") and not hasattr(method, "resolve_boundary_conditions")
    with pytest.raises(B.UnboundTagError):
        H.diffusion_mesh(H.chart_mesh("slab", left=0, right=0))

r"""B1 regression floor — the affine form's two factors are POPULATED.

Campaign phase B1 (``.claude/plans/boundary_machinery_review.md``) minted
:class:`~orpheus.geometry.boundary.BoundaryGeometryMap` /
:class:`~orpheus.geometry.boundary.BoundaryResponseKernel` and populated them on
all seven production laws. Before it, both were ``-> Any`` properties defaulting
to ``None`` that **nothing populated and nothing read**, while five production
sites answered the questions they exist to answer by comparing ``law.kind``
strings.

This file is the floor that keeps them populated. The ABC default is still
``None`` (a stub law inherits it — pinned in
:mod:`tests.geometry.test_boundary_trace_law`), so nothing structurally stops a
future law from shipping unpopulated; these tests are what would catch it.

The bit-identity leg is the phase's stated acceptance criterion: the minted
``response_kernel`` must be the **same number** the realizer already reaches
for, not a re-derivation of it. A new number here would mean B1 changed physics
while claiming to be a pure addition. (There were TWO such legs; B2.2 retired
the diffusion one when the collapse made it tautological — see the retirement
note in the body, which is the worked example of that decay.)

Tagged ``@pytest.mark.foundation``: these assert a structural contract, not a
discretization claim, so they carry no ``verifies(...)`` equation label (the
verifies⊥level doctrine).
"""

from __future__ import annotations

import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryGeometryMap,
    BoundaryResponseKernel,
    BoundaryTraceLaw,
    HemisphericalAverage,
    IdentityMap,
    NullMap,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    ScalarResponse,
    SpatialWrap,
    SpecularMirror,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)


#: ``(law, expected geometry type, expected response scalar, test id)`` — one
#: entry per production law. The parametrization below AND the
#: registry-coverage check both derive from this single list, so they cannot
#: drift apart.
#:
#: ``_stub_for_test`` is deliberately absent: it is a test fixture, and the
#: registry-hygiene fixture in ``test_boundary_trace_law.py`` evicts it outside
#: that module anyway.
_LAW_SPECS: list[tuple[BoundaryTraceLaw, type, float, str]] = [
    (VacuumInflow(), NullMap, 0.0, "vacuum"),
    (ReflectiveBoundary(), SpecularMirror, 1.0, "reflective-a1"),
    (ReflectiveBoundary(axis="y", albedo=0.7), SpecularMirror, 0.7,
     "reflective-partial"),
    (WhiteBoundary(axis="x", albedo=0.4), HemisphericalAverage, 0.4, "white"),
    (AlbedoBoundary(albedo=0.25), IdentityMap, 0.25, "albedo"),
    (PeriodicBoundary(), SpatialWrap, 1.0, "periodic"),
    (PrescribedInflow(), NullMap, 0.0, "prescribed_inflow"),
    (ZeroFluxBoundary(), IdentityMap, -1.0, "zero_flux"),
]

PRODUCTION_LAWS = [
    pytest.param(law, geom, alpha, id=tid)
    for law, geom, alpha, tid in _LAW_SPECS
]


@pytest.mark.foundation
@pytest.mark.parametrize("law,geom_cls,alpha", PRODUCTION_LAWS)
def test_every_production_law_states_both_factors(law, geom_cls, alpha) -> None:
    """No production law may report ``None`` for either factor.

    This is the whole of B1 as a single assertion. It reds if a new law ships
    without its spec, or if someone reverts a population.
    """
    if law.geometry_map is None or law.response_kernel is None:
        pytest.fail(
            f"{type(law).__name__} reports geometry_map="
            f"{law.geometry_map!r}, response_kernel={law.response_kernel!r} — "
            f"an unpopulated factor is the pre-B1 state that forced production "
            f"to dispatch on `law.kind` strings instead."
        )
    assert isinstance(law.geometry_map, geom_cls)
    assert isinstance(law.response_kernel, ScalarResponse)


@pytest.mark.foundation
@pytest.mark.parametrize("law,geom_cls,alpha", PRODUCTION_LAWS)
def test_factors_satisfy_their_protocols(law, geom_cls, alpha) -> None:
    """Structural conformance to the two Protocols, at runtime."""
    assert isinstance(law.geometry_map, BoundaryGeometryMap)
    assert isinstance(law.response_kernel, BoundaryResponseKernel)


@pytest.mark.foundation
@pytest.mark.parametrize("law,geom_cls,alpha", PRODUCTION_LAWS)
def test_response_scalar_is_the_declared_amplitude(law, geom_cls, alpha) -> None:
    """``response_kernel.scalar`` is exactly the law's amplitude."""
    scalar = law.response_kernel.scalar
    if scalar != alpha:
        pytest.fail(
            f"{type(law).__name__}.response_kernel.scalar = {scalar!r}, "
            f"expected exactly {alpha!r}. B1 is a pure addition — the minted "
            f"factor must NAME the existing number, never re-derive it."
        )


@pytest.mark.foundation
@pytest.mark.parametrize("law,geom_cls,alpha", PRODUCTION_LAWS)
def test_response_is_bit_identical_to_the_sn_realizer_float(
    law, geom_cls, alpha,
) -> None:
    r"""The SN arm multiplies by ``float(law.albedo)``; we must equal it.

    ``sn/boundary/realizer.py`` builds ``float(law.albedo) * base`` for every
    :math:`\alpha \ne 1` law. If the minted scalar diverged from that, B1 would
    have changed the realized operator while claiming to touch nothing.
    """
    if not hasattr(law, "albedo"):
        pytest.skip(
            f"{type(law).__name__} carries no `albedo` field — its response is "
            f"structural (0 for the rank-0 laws, 1 for periodic), so there is "
            f"no realizer float to compare against."
        )
    assert law.response_kernel.scalar == float(law.albedo)


# ── RETIRED at B2.2, and why ──────────────────────────────────────────
#
# ``test_response_is_bit_identical_to_the_diffusion_partial_current_albedo``
# lived here. It asserted
#
#     law.response_kernel.scalar == DiffusionBoundaryRealizer\
#                                       ._partial_current_albedo(law)
#
# which was real evidence in B1: the realizer was a five-arm ``isinstance``
# ladder returning ``0.0`` / ``-1.0`` / ``float(law.albedo)``, and the equality
# proved the minted factor was a rename-and-lift of the number that ladder
# already reached for, not a re-derivation.
#
# B2.2 collapsed the ladder onto ``law.response_kernel.scalar``. The assertion
# is now ``x == x`` — it executes, it is green, and it can never red: the same
# fires-but-cannot-fail family as a tautological companion guard
# (``vv-principles`` Mode 8), arrived at not by bad authorship but by DECAY,
# because production moved underneath a correctly-written test.
#
# Its coverage is MIGRATED, not dropped (coding-standards: retirement means
# test migration). ``tests/geometry/test_boundary_factor_consumers.py``
# reproduces the retired ladder verbatim and compares it against the live
# realizer law by law — the same claim, with the teeth restored, plus a
# mutation leg proving 𝒜 tracks the factor.
#
# The SN sibling ABOVE is NOT retired: the SN realizer still spells its own
# ``float(law.albedo) * base``, so that comparison still has two sides. It
# becomes tautological at phase B4, and should be retired the same way then.


@pytest.mark.foundation
def test_every_registered_law_is_covered_by_this_file() -> None:
    """The parametrization must not silently fall behind the registry.

    Without this, adding an eighth law would leave it untested here while every
    test above still passed — the coverage would look complete and be stale.
    """
    registered = {
        k for k in BoundaryTraceLaw.registry if not k.startswith("_")
    }
    # ``kind`` is the registry key for every law EXCEPT a partially-reflecting
    # ReflectiveBoundary, which reports "partial" (see the B2 warning below);
    # compare on the registered class instead so this check is about coverage,
    # not about that wrinkle.
    covered = {type(law).key for law, _, _, _ in _LAW_SPECS}
    missing = registered - covered
    if missing:
        pytest.fail(
            f"law(s) {sorted(missing)} are registered but absent from "
            f"PRODUCTION_LAWS — add them, with their factors."
        )


@pytest.mark.foundation
@pytest.mark.parametrize("law,geom_cls,alpha", PRODUCTION_LAWS)
def test_factors_are_frozen(law, geom_cls, alpha) -> None:
    """A spec is a value; mutating one in place must be impossible."""
    import dataclasses

    for factor in (law.geometry_map, law.response_kernel):
        fields = dataclasses.fields(factor)
        if not fields:
            # IdentityMap / NullMap are parameterless — there is nothing to
            # mutate, so frozen-ness is vacuous for them rather than untested.
            continue
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(factor, fields[0].name, None)


@pytest.mark.foundation
def test_specular_mirror_is_the_only_ordinate_permuting_geometry() -> None:
    r"""``permutes_ordinates`` is the structural form of ``== "reflective"``.

    ``sweep_schedule.py`` currently selects reflective faces by string compare;
    phase B2 repoints it here. This pins the answer set so that repointing is a
    behaviour-preserving change rather than a redefinition.

    White answers ``False`` deliberately: it couples every outgoing ordinate to
    every incoming one (rank-1 in angle), which is an all-to-all coupling, not
    a relabeling — and is why the production Gauss-Seidel schedule excludes
    white from its octant split.

    .. note::

       **A correction to this file's own first draft.** It warned that
       repointing ``sweep_schedule.py`` here would NOT be behaviour-preserving,
       reasoning that ``ReflectiveBoundary(albedo=0.7).kind`` is ``"partial"``
       so the string compare must miss partial reflectors. **Measured, that is
       wrong** — production never reads ``law.kind``:

       * ``SNMesh.realize_boundary_law`` returns
         ``_BoundBoundaryOperator(realized, kind=law.key)``
         (``sn/mesh/augmented_mesh.py:435``) — it stores the **registry key**,
       * and ``ReflectiveBoundary.key`` is ``"reflective"`` for **every**
         albedo, including 0.7.

       So ``bc[face] == "reflective"`` already matches partial reflectors, and
       :attr:`permutes_ordinates` is ``True`` for both — **they agree**, and
       repointing is behaviour-preserving. The legs below pin that agreement so
       B2 can rely on it.

       The wrinkle that IS real: ``law.kind`` and the wrapper's stored ``key``
       diverge for a partial reflector, and only the latter reaches production
       — because the wrapper discards the law entirely.
    """
    permuting_geoms = {
        type(law.geometry_map)
        for law, _, _, _ in _LAW_SPECS
        if law.geometry_map.permutes_ordinates
    }
    assert permuting_geoms == {SpecularMirror}, (
        f"expected exactly SpecularMirror to permute ordinates; got "
        f"{sorted(c.__name__ for c in permuting_geoms)}. A change here "
        f"silently alters which faces the sweep schedule treats as reflective."
    )

    # The agreement B2 relies on, pinned. A partial reflector permutes AND
    # carries the "reflective" registry key that production actually compares
    # against, so the structural predicate and the string compare give the same
    # answer for every reachable face.
    partial = ReflectiveBoundary(axis="y", albedo=0.7)
    assert partial.geometry_map.permutes_ordinates is True
    assert type(partial).key == "reflective"   # what the wrapper stores
    assert partial.kind == "partial"           # what the LAW says — divergent,
    #                                            but production never sees it

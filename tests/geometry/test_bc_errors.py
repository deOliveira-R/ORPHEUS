r"""Foundation tests for the 8 named boundary-condition errors.

This file pins BOTH halves of each typed error's contract:

**(1) The TYPE contract** — on a hand-built instance:

* ``isinstance(err, BoundaryError)`` and the inherited
  ``ValueError`` MRO (so ``except ValueError`` consumers continue
  to catch boundary violations).
* ``err.law`` carries the offending-law name passed at construction.
* ``str(err)`` yields the message passed at construction.

**(2) The PRODUCTION contract** — the guard that actually raises it:

each test drives a **production entry point** (a law invariant, or
``SNBoundaryRealizer.realize`` / ``DiffusionBoundaryRealizer.realize``)
into the violation and asserts the escaping exception is of the right
type AND carries the right ``law=`` tag.

Why (2) exists — the B0.3 repair
--------------------------------

Until 2026-07-30 every leg here read

.. code-block:: python

    err = SomeBoundaryError("msg", law="x")
    with pytest.raises(SomeBoundaryError):
        raise err            # <- err was built two lines above AS a SomeBoundaryError

That is a **tautology** (``vv-principles`` Mode 8, the SELF-SATISFIED
``pytest.raises`` class): it verifies that Python raises what you told
it to, and pins nothing about production. Measured in the boundary
review: **zero of 14 deliberate guard-disabling mutations reddened this
file**, while every one of them reddened the real negative tests in
``test_bc_universal_invariants.py`` / ``test_sn_boundary_realizer.py``.
The file read as "the error paths are heavily tested" while contributing
zero mutation detection.

The acceptance test for the repair is stated as a mutation: **disable
the production guard named in a test's docstring and THIS file must
red.** Each test names its guard at ``file:line`` so that check is
mechanical.

One error has no production raiser
----------------------------------

**[G]** ``IncomingOutgoingTraceClassificationError`` (ERR-040) has
**zero** ``raise`` sites in ``orpheus/`` — its documented trigger,
``BoundaryTraceLaw.assert_inflow_outflow_classification``
(``_base.py:251``), is a no-op ABC default that **no** concrete law
overrides. It is therefore type-contract-only below, deliberately and
with the absence recorded, rather than wearing a fake production leg.

Tagged ``@pytest.mark.foundation`` (software invariant on the
type + its guard wiring, not an L0 claim against a theory page) per
:mod:`tests._harness` conventions.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    ConstantInflowSource,
    IncomingOutgoingTraceClassificationError,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectionNotInvolutiveError,
    ReflectiveBoundary,
    SubmarkovViolationError,
    VacuumAppliedToOutgoingTraceError,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace


class _TableQuad:
    """Quadrature stand-in with an INJECTABLE reflection table.

    Only the surface ``ReflectiveBoundary.assert_is_involutive`` reads
    (``N`` + ``reflection_index``); every other datum is irrelevant to
    that guard. Used ONLY where the realize-time route cannot isolate
    the invariant under test (see
    ``test_reflection_not_involutive_error``).
    """

    def __init__(self, n: int, table: np.ndarray) -> None:
        self.N = n
        self._table = table

    def reflection_index(self, axis: str) -> np.ndarray:  # noqa: ARG002
        return self._table


# ─────────────────────────────────────────────────────────────────────
# Type contract + production guard, per error class
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_boundary_error_base_class_contract() -> None:
    """:class:`BoundaryError` is a :class:`ValueError` subclass, and
    both realizers raise it bare on an un-dispatchable law.

    Existing ``except ValueError`` consumers in the boundary layer
    must continue to catch :class:`BoundaryError` and every typed
    subclass -- this is the backward-compatibility hook.

    PRODUCTION GUARDS: ``sn/boundary/realizer.py:171`` (SN refuses
    ``ZeroFluxBoundary`` -- a negative angular inflow has no transport
    realization), ``diffusion/boundary_realizer.py:219`` (diffusion
    refuses ``PeriodicBoundary`` -- a trace-block wrap is not a
    per-face albedo) and ``diffusion/boundary_realizer.py:229``
    (diffusion refuses ``PrescribedInflow`` -- the rank-0 AFFINE law
    realizes as a source, not a linear operator).
    """
    err = BoundaryError("msg", law="my_law")
    assert isinstance(err, ValueError)
    assert isinstance(err, BoundaryError)
    assert err.law == "my_law"
    assert str(err) == "msg"

    quad = Quadrature.lebedev(17)
    with pytest.raises(BoundaryError) as sn_exc:
        SNBoundaryRealizer().realize(
            ZeroFluxBoundary(), SNMethodSpace.minimal(quad),
        )
    assert sn_exc.value.law == "zero_flux"
    assert isinstance(sn_exc.value, ValueError)

    with pytest.raises(BoundaryError) as periodic_exc:
        DiffusionBoundaryRealizer().realize(
            PeriodicBoundary(), DiffusionMethodSpace.minimal(),
        )
    assert periodic_exc.value.law == "periodic"

    with pytest.raises(BoundaryError) as prescribed_exc:
        DiffusionBoundaryRealizer().realize(
            PrescribedInflow(source=ConstantInflowSource(value=1.0)),
            DiffusionMethodSpace.minimal(),
        )
    assert prescribed_exc.value.law == "prescribed_inflow"

    # Positive control: a law BOTH realizers accept must NOT raise, so
    # the three reds above come from the refusal branches and not from
    # an always-raising realizer.
    SNBoundaryRealizer().realize(
        ReflectiveBoundary(axis="x"), SNMethodSpace.minimal(quad),
    )
    DiffusionBoundaryRealizer().realize(
        ReflectiveBoundary(axis="x"), DiffusionMethodSpace.minimal(),
    )


@pytest.mark.foundation
def test_incoming_outgoing_trace_classification_error() -> None:
    """TYPE CONTRACT ONLY -- this error has **no** production raiser.

    **[G]** ``grep -rn "raise IncomingOutgoingTraceClassificationError"
    orpheus/`` returns zero hits (2026-07-30). Its documented trigger,
    :meth:`BoundaryTraceLaw.assert_inflow_outflow_classification`
    (``orpheus/geometry/boundary/_base.py:251``), is a no-op ABC
    default that **no** concrete law overrides, so there is no
    production entry point to point a ``pytest.raises`` leg at.

    The honest gate is therefore the type contract alone. Deliberately
    NOT given a synthetic raise leg: a ``with pytest.raises(X): raise
    err`` here would be exactly the tautology this file was repaired to
    remove. If a law ever overrides the classification invariant, add
    the production leg then.
    """
    err = IncomingOutgoingTraceClassificationError(
        "tangential ordinate at face 'left'", law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    assert str(err) == "tangential ordinate at face 'left'"


@pytest.mark.foundation
def test_vacuum_applied_to_outgoing_trace_error() -> None:
    """ERR-041: vacuum realized against an OUTGOING trace annotation.

    PRODUCTION GUARD: ``sn/boundary/realizer.py:237`` -- the realizer
    cross-checks the claimed ``inflow_indices`` against the signed
    projection the FACE NAME implies, and refuses any claimed inflow
    ordinate that is outgoing. Here ``xmax``'s outgoing set
    (:math:`\\mu_x > 0`) is handed in as its inflow set.
    """
    err = VacuumAppliedToOutgoingTraceError(
        "vacuum BC on Gamma_+", law="vacuum",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "vacuum"
    assert str(err) == "vacuum BC on Gamma_+"

    quad = Quadrature.gauss_legendre(n_ordinates=8)
    swapped = SNMethodSpace(
        quadrature=quad, face="xmax",
        inflow_indices=np.flatnonzero(quad.mu_x > 0),
    )
    with pytest.raises(VacuumAppliedToOutgoingTraceError) as exc:
        SNBoundaryRealizer().realize(VacuumInflow(), swapped)
    assert exc.value.law == "vacuum"


@pytest.mark.foundation
def test_boundary_geometry_map_not_measure_preserving_error(monkeypatch) -> None:
    """ERR-042: a reflection table that does not preserve
    :math:`w\\cdot|\\mu|`.

    PRODUCTION GUARD:
    ``ReflectiveBoundary.assert_geometry_map_measure_preserving``
    (``geometry/boundary/reflective.py:212``), fired through
    ``SNBoundaryRealizer.realize`` -> ``assert_realizable``. The
    neighbour-pair table on GL-8 is a perfect involution that pairs
    ordinates of UNEQUAL ``w·|mu|`` -- only the measure check sees it.
    """
    err = BoundaryGeometryMapNotMeasurePreservingError(
        "reflection table inconsistent with quadrature weights",
        law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    assert str(err) == (
        "reflection table inconsistent with quadrature weights"
    )

    quad = Quadrature.gauss_legendre(n_ordinates=8)
    monkeypatch.setitem(
        quad.reflection_partners, 0,
        np.array([1, 0, 3, 2, 5, 4, 7, 6]),
    )
    with pytest.raises(
        BoundaryGeometryMapNotMeasurePreservingError,
    ) as exc:
        SNBoundaryRealizer().realize(
            ReflectiveBoundary(axis="x"), SNMethodSpace.minimal(quad),
        )
    assert exc.value.law == "reflective"


@pytest.mark.foundation
def test_boundary_response_not_positive_error() -> None:
    """ERR-043: a negative response amplitude.

    PRODUCTION GUARDS:
    ``WhiteBoundary.assert_response_positive_if_declared``
    (``geometry/boundary/white.py:101``) and
    ``AlbedoBoundary.assert_response_positive_if_declared``
    (``geometry/boundary/albedo.py:78``). Both are fired by the
    laws' own ``assert_realizable`` at the realizer seam; the direct
    invariant call is used here so each guard is attributable
    one-to-one.
    """
    err = BoundaryResponseNotPositiveError(
        "white kernel produced negative output",
        law="white",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "white"

    with pytest.raises(BoundaryResponseNotPositiveError) as white_exc:
        WhiteBoundary(
            axis="x", outward_sign=+1, albedo=-0.1,
        ).assert_response_positive_if_declared()
    assert white_exc.value.law == "white"

    with pytest.raises(BoundaryResponseNotPositiveError) as albedo_exc:
        AlbedoBoundary(albedo=-0.1).assert_response_positive_if_declared()
    assert albedo_exc.value.law == "albedo"


@pytest.mark.foundation
def test_reflection_not_involutive_error() -> None:
    """ERR-044: ``perm[perm] != arange(N)``.

    PRODUCTION GUARD: ``ReflectiveBoundary.assert_is_involutive``
    (``geometry/boundary/reflective.py:155``).

    Driven by the DIRECT invariant call rather than through
    ``realize``: a rotation table also breaks the ``w·|mu|`` measure,
    and ``assert_realizable`` fires the measure check FIRST
    (``reflective.py:290`` calls ``super()`` before the involution),
    so the realize-time route would attribute the red to ERR-042. The
    guard under test is production code either way -- only the
    quadrature carrying the bad table is scaffolding.
    """
    err = ReflectionNotInvolutiveError(
        "perm[perm] != arange(N)",
        law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    assert str(err) == "perm[perm] != arange(N)"

    n = 8
    rotation = _TableQuad(n, np.roll(np.arange(n), 1))
    with pytest.raises(ReflectionNotInvolutiveError) as exc:
        ReflectiveBoundary(axis="x").assert_is_involutive(rotation)  # type: ignore[arg-type]
    assert exc.value.law == "reflective"

    # Positive control: a true involution on the SAME shim passes, so
    # the red above is attributable to the table, not the stand-in.
    involution = _TableQuad(n, np.array([1, 0, 3, 2, 5, 4, 7, 6]))
    ReflectiveBoundary(axis="x").assert_is_involutive(involution)  # type: ignore[arg-type]


@pytest.mark.foundation
def test_reflection_did_not_map_inflow_to_outflow_error(monkeypatch) -> None:
    """ERR-045: a reflection that maps an inflow ordinate to itself.

    PRODUCTION GUARD:
    ``ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow``
    (``geometry/boundary/reflective.py:264``), fired through
    ``SNBoundaryRealizer.realize`` -> ``assert_realizable``. The
    identity table is the canonical mutant: trivially involutive AND
    trivially measure-preserving, so it sails past the two sibling
    invariants and only this one sees it.
    """
    err = ReflectionDidNotMapInflowToOutflowError(
        "inflow ordinate maps to inflow under reflection",
        law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    assert str(err) == (
        "inflow ordinate maps to inflow under reflection"
    )

    quad = Quadrature.gauss_legendre(n_ordinates=8)
    monkeypatch.setitem(
        quad.reflection_partners, 0, np.arange(quad.N),
    )
    with pytest.raises(ReflectionDidNotMapInflowToOutflowError) as exc:
        SNBoundaryRealizer().realize(
            ReflectiveBoundary(axis="x"), SNMethodSpace.minimal(quad),
        )
    assert exc.value.law == "reflective"


@pytest.mark.foundation
def test_submarkov_violation_error() -> None:
    """ERR-046: :math:`\\alpha > 1` on a sub-Markov BC.

    PRODUCTION GUARDS: ``AlbedoBoundary.assert_submarkov``
    (``geometry/boundary/albedo.py:93``) and
    ``WhiteBoundary.assert_submarkov``
    (``geometry/boundary/white.py:117``).
    """
    err = SubmarkovViolationError(
        "albedo=1.2 violates row-sum <= 1",
        law="albedo",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "albedo"

    with pytest.raises(SubmarkovViolationError) as albedo_exc:
        AlbedoBoundary(albedo=1.2).assert_submarkov()
    assert albedo_exc.value.law == "albedo"

    with pytest.raises(SubmarkovViolationError) as white_exc:
        WhiteBoundary(
            axis="x", outward_sign=+1, albedo=1.2,
        ).assert_submarkov()
    assert white_exc.value.law == "white"

    # Positive control (vv-principles anti-pattern #11): the same
    # guards MUST NOT raise on a legal amplitude, so the reds above
    # come from the bound, not from an always-raising body.
    AlbedoBoundary(albedo=1.0).assert_submarkov()
    WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0).assert_submarkov()


@pytest.mark.foundation
def test_boundary_source_not_on_incoming_trace_error() -> None:
    """ERR-047: an affine source with support off :math:`\\Gamma_-`.

    PRODUCTION GUARD:
    ``BoundaryTraceLaw.assert_source_lives_on_incoming_trace``
    (``geometry/boundary/_base.py:333``), fired both directly and
    through ``SNBoundaryRealizer.realize``: a nonzero
    ``ConstantInflowSource`` on a method space with NO inflow indices
    has no mask to confine it to the incoming trace.
    """
    err = BoundarySourceNotOnIncomingTraceError(
        "q has nonzero entries on outflow ordinates",
        law="prescribed_inflow",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "prescribed_inflow"
    assert str(err) == "q has nonzero entries on outflow ordinates"

    quad = Quadrature.lebedev(17)
    law = PrescribedInflow(source=ConstantInflowSource(value=2.0))
    with pytest.raises(BoundarySourceNotOnIncomingTraceError) as exc:
        SNBoundaryRealizer().realize(law, SNMethodSpace.minimal(quad))
    # B0.3 FINDING (measured): every other guard tags ``law=`` with the
    # lowercase REGISTRY KEY -- "vacuum", "reflective", "white",
    # "albedo", "zero_flux", "periodic" (all asserted verbatim above).
    # This one alone tags ``law=type(self).__name__``
    # (``geometry/boundary/_base.py:342``), so the escaping tag reads
    # "PrescribedInflow", not "prescribed_inflow". Asserted
    # case/underscore-insensitively so the leg pins that the tag
    # IDENTIFIES the law -- dropping the ``law=`` kwarg still reds it --
    # WITHOUT calcifying the drift the B0.3 report asks to fix at the
    # production site.
    assert exc.value.law.lower().replace("_", "") == "prescribedinflow"

    # Positive control: with the face's inflow indices supplied the
    # same guard certifies silently.
    law.assert_source_lives_on_incoming_trace(
        quad, np.flatnonzero(quad.mu_x < 0),
    )


# ─────────────────────────────────────────────────────────────────────
# Default-argument contract
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_law_field_defaults_to_empty_string() -> None:
    """``law`` is an optional kwarg; default is the empty string.

    Some construction sites (e.g. invariant checks raised before a
    concrete law has been identified) won't carry a law name. The
    error must still be raise-able.
    """
    err = BoundaryError("anonymous violation")
    assert err.law == ""
    assert str(err) == "anonymous violation"


@pytest.mark.foundation
def test_typed_subclasses_inherit_law_field_default() -> None:
    """Each typed subclass inherits ``BoundaryError``'s ``law=``
    default, so a caller can choose to omit it without breaking the
    constructor."""
    err = BoundarySourceNotOnIncomingTraceError("anonymous")
    assert err.law == ""
    assert isinstance(err, BoundaryError)

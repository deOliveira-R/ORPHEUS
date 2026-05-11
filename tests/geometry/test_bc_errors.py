r"""Foundation tests for the 8 named boundary-condition errors.

This file pins the TYPE contract for the
:mod:`orpheus.geometry.boundary._errors` module. Each test
instantiates one error class, verifies:

* ``isinstance(err, BoundaryError)`` and the inherited
  ``ValueError`` MRO (so ``except ValueError`` consumers continue
  to catch boundary violations).
* ``err.law`` carries the offending-law name passed at construction.
* ``str(err)`` yields the message passed at construction.
* The error is ``raise``-able and ``pytest.raises`` selects it
  correctly.

Tagged ``@pytest.mark.foundation`` (software invariant on the
type, not an L0 claim against a theory page) per
:mod:`tests._harness` conventions.

The ``@pytest.mark.catches("ERR-NNN")`` decorators that wire each
error class to its catalogue entry are NOT added here -- they
ship in Wave 7 alongside the concrete BCs that fire the errors
through invariant overrides. Wave 3's job is to pin the TYPES.
"""

from __future__ import annotations

import pytest

from orpheus.geometry.boundary import (
    BoundaryError,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    IncomingOutgoingTraceClassificationError,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectionNotInvolutiveError,
    SubmarkovViolationError,
    VacuumAppliedToOutgoingTraceError,
)


# ─────────────────────────────────────────────────────────────────────
# Type contract per error class
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_boundary_error_base_class_contract() -> None:
    """:class:`BoundaryError` is a :class:`ValueError` subclass.

    Existing ``except ValueError`` consumers in the boundary layer
    must continue to catch :class:`BoundaryError` and every typed
    subclass -- this is the backward-compatibility hook.
    """
    err = BoundaryError("msg", law="my_law")
    assert isinstance(err, ValueError)
    assert isinstance(err, BoundaryError)
    assert err.law == "my_law"
    assert str(err) == "msg"
    with pytest.raises(BoundaryError):
        raise err


@pytest.mark.foundation
def test_incoming_outgoing_trace_classification_error() -> None:
    err = IncomingOutgoingTraceClassificationError(
        "tangential ordinate at face 'left'", law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    assert str(err) == "tangential ordinate at face 'left'"
    with pytest.raises(IncomingOutgoingTraceClassificationError):
        raise err


@pytest.mark.foundation
def test_vacuum_applied_to_outgoing_trace_error() -> None:
    err = VacuumAppliedToOutgoingTraceError(
        "vacuum BC on Gamma_+", law="vacuum",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "vacuum"
    assert str(err) == "vacuum BC on Gamma_+"
    with pytest.raises(VacuumAppliedToOutgoingTraceError):
        raise err


@pytest.mark.foundation
def test_boundary_geometry_map_not_measure_preserving_error() -> None:
    err = BoundaryGeometryMapNotMeasurePreservingError(
        "reflection table inconsistent with quadrature weights",
        law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError):
        raise err


@pytest.mark.foundation
def test_boundary_response_not_positive_error() -> None:
    err = BoundaryResponseNotPositiveError(
        "white kernel produced negative output",
        law="white",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "white"
    with pytest.raises(BoundaryResponseNotPositiveError):
        raise err


@pytest.mark.foundation
def test_reflection_not_involutive_error() -> None:
    err = ReflectionNotInvolutiveError(
        "perm[perm] != arange(N)",
        law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    with pytest.raises(ReflectionNotInvolutiveError):
        raise err


@pytest.mark.foundation
def test_reflection_did_not_map_inflow_to_outflow_error() -> None:
    err = ReflectionDidNotMapInflowToOutflowError(
        "inflow ordinate maps to inflow under reflection",
        law="reflective",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "reflective"
    with pytest.raises(ReflectionDidNotMapInflowToOutflowError):
        raise err


@pytest.mark.foundation
def test_submarkov_violation_error() -> None:
    err = SubmarkovViolationError(
        "albedo=1.2 violates row-sum <= 1",
        law="albedo",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "albedo"
    with pytest.raises(SubmarkovViolationError):
        raise err


@pytest.mark.foundation
def test_boundary_source_not_on_incoming_trace_error() -> None:
    err = BoundarySourceNotOnIncomingTraceError(
        "q has nonzero entries on outflow ordinates",
        law="prescribed_inflow",
    )
    assert isinstance(err, BoundaryError)
    assert isinstance(err, ValueError)
    assert err.law == "prescribed_inflow"
    assert str(err) == "q has nonzero entries on outflow ordinates"
    with pytest.raises(BoundarySourceNotOnIncomingTraceError):
        raise err


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
    default, so a Wave-7 caller can choose to omit it without
    breaking the constructor."""
    err = BoundarySourceNotOnIncomingTraceError("anonymous")
    assert err.law == ""
    assert isinstance(err, BoundaryError)

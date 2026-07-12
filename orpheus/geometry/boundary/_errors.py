r"""Named-error catalog for boundary-condition violations.

Per Grand Report v3 §26A.4 lines 5072-5095, each typed error replaces
a generic :class:`ValueError` raised by a :class:`BoundaryTraceLaw`
invariant or realisation step. Every error class is mapped to an
ERR-NNN entry in ``.claude/skills/vv-principles/error_catalog.md``
and to a ``@pytest.mark.catches("ERR-NNN")`` decorator on the test
that proves the error fires under the right conditions.

The 8 typed errors below replace generic :class:`ValueError` raises
in the boundary layer. Wave 3 shipped the TYPES; Wave 7 shipped the
concrete BCs firing the response/involution/sub-Markov family; #52
(2026-07-12) wired the remaining four as PRODUCTION guards — the
measure-preservation direct check (ERR-042), the inflow→outflow
table check (ERR-045), the universal source-trace certification
(ERR-047) all fire through
:meth:`~orpheus.geometry.boundary._base.BoundaryTraceLaw.assert_realizable`
at every realize, and the vacuum trace-orientation guard (ERR-041)
fires at the SN realizer's vacuum arm.

References
----------

* Grand Report v3 §16A.12 (universal invariants) and §27.6
  (named-error catalog discipline).
* ``.claude/plans/transient-giggling-cake.md`` — Wave 3 brief.
"""

from __future__ import annotations

__all__ = [
    "BoundaryError",
    "IncomingOutgoingTraceClassificationError",
    "VacuumAppliedToOutgoingTraceError",
    "BoundaryGeometryMapNotMeasurePreservingError",
    "BoundaryResponseNotPositiveError",
    "ReflectionNotInvolutiveError",
    "ReflectionDidNotMapInflowToOutflowError",
    "SubmarkovViolationError",
    "BoundarySourceNotOnIncomingTraceError",
]


class BoundaryError(ValueError):
    """Base class for boundary-condition violation errors.

    Carries a ``law`` field naming which :class:`BoundaryTraceLaw`
    subclass raised the violation, so consumers can pattern-match
    on the error class and the offending law in one step. Extends
    :class:`ValueError` for backward compatibility with existing
    ``except ValueError`` consumers in the boundary layer.

    Parameters
    ----------
    message : str
        Human-readable message describing the violation.
    law : str, optional
        Name of the :class:`BoundaryTraceLaw` subclass that raised
        the violation. Default is the empty string for cases where
        the law is unknown or doesn't apply.
    """

    def __init__(self, message: str, law: str = "") -> None:
        super().__init__(message)
        self.law = law


class IncomingOutgoingTraceClassificationError(BoundaryError):
    r"""An ordinate could not be classified as inflow or outflow.

    Raised when ``assert_inflow_outflow_classification`` finds an
    ordinate with :math:`|\Omega \cdot \hat n| \leq \epsilon` on a
    face where a strict inflow/outflow partition was required (no
    tangential ordinates allowed by the law's contract).
    """


class VacuumAppliedToOutgoingTraceError(BoundaryError):
    r"""A vacuum BC was constructed against an outgoing trace.

    Vacuum sets :math:`\gamma_- \psi = 0` on the INCOMING trace
    only. Applying it to an outgoing trace is geometrically
    meaningless and typically indicates a mesh/face annotation
    bug.
    """


class BoundaryGeometryMapNotMeasurePreservingError(BoundaryError):
    r"""The geometric map :math:`G_\alpha` is not measure-preserving.

    For permutation-based maps (specular reflection) the
    direction-cosine measure :math:`w(\Omega)\,|\Omega\cdot n|`
    must be invariant under the map. Violation indicates a wrong
    reflection-index table or a quadrature whose ``mu_n`` is
    inconsistent with its weights.
    """


class BoundaryResponseNotPositiveError(BoundaryError):
    r"""The response kernel :math:`R` produces negative values.

    For physically meaningful BCs (albedo, white) the response is
    non-negative on the inflow trace. Negative outputs signal a
    sign error in the kernel construction.
    """


class ReflectionNotInvolutiveError(BoundaryError):
    r"""A reflection permutation is not an involution: :math:`\pi \circ \pi \neq \mathrm{id}`.

    For a clean axis reflection, applying the reflection twice
    must return every ordinate to itself. A non-involutive
    reflection table is typically a bug in the reflection-index
    construction (wrong axis, wrong ordering, or partial
    permutation).
    """


class ReflectionDidNotMapInflowToOutflowError(BoundaryError):
    r"""A reflection mapped an inflow ordinate to itself (not outflow).

    For a clean axis reflection, every inflow ordinate at the face
    maps to an outflow ordinate (and vice versa). Violation
    indicates a wrong reflection-index table or a non-axis-aligned
    reflection that requires a different BC type (``SymmetryBoundary``
    in the future) rather than ``ReflectiveBoundary``.
    """


class SubmarkovViolationError(BoundaryError):
    r"""A sub-Markov BC (albedo, white) has :math:`\alpha > 1`.

    A sub-Markov kernel satisfies :math:`\int R\,\mathrm{d}y \leq 1`
    on each row -- equivalently, the total reflected current does
    not exceed the incident current. Albedo :math:`> 1` violates
    this physically (would imply a source on the boundary surface).
    """


class BoundarySourceNotOnIncomingTraceError(BoundaryError):
    r"""A boundary source :math:`q` has nonzero entries on the outflow trace.

    The affine BC form
    :math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q` requires
    :math:`q \in \Gamma_-` (incoming trace only). Outflow-trace
    sources are geometrically meaningless and indicate a wrong
    source-shape contract.
    """

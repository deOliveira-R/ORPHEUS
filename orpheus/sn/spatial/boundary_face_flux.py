r"""Boundary face-flux strategy contract for the curvilinear FD operator.

Why this abstraction exists
===========================

The symmetric-closure transport operator
:func:`~orpheus.sn.operator.transport_operator_matvec_spherical` /
:func:`~orpheus.sn.operator.transport_operator_matvec_cylindrical` reads
the radial face flux :math:`\psi_{i \pm 1/2}` at every cell.  In the
**interior** the symmetric arithmetic average

.. math::

   \psi_{i+1/2} \;\approx\; \tfrac{1}{2}\,
       \bigl(\psi_i + \psi_{i+1}\bigr)

is :math:`\mathcal{O}(h^2)` accurate.  At the **boundary cell**
:math:`i = N - 1` for an outgoing ordinate (:math:`\mu > 0`), there is
no neighbour cell at :math:`i = N`; the historical operator used the
cell-centre value directly,

.. math::

   \psi_{N-1/2}^{\rm out} \;\approx\; \psi_{N-1},

which is only :math:`\mathcal{O}(h)` accurate on a non-constant
solution and was Issue #168 Defect 1: the leading boundary truncation
that limits the curvilinear MMS convergence rate.  Phase A of Option C
replaces this cell-centre closure with the **DD diamond extrapolation**

.. math::
   :label: dd-boundary-extrapolation

   \psi_{N-1/2}^{\rm out} \;=\; \tfrac{3}{2}\,\psi_{N-1}
       \;-\; \tfrac{1}{2}\,\psi_{N-2}

— a one-sided second-order extrapolation through two consecutive cell
centres (equivalently, the DD relation
:math:`\psi^{\rm face}_{i+1/2} = 2\,\bar\psi_i -
\psi^{\rm face}_{i-1/2}` solved for the outer face given the
interior face value :math:`\tfrac{1}{2}(\psi_{N-2} + \psi_{N-1})`).

Lifting the boundary closure into a strategy class — analogous to the
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` Protocol —  makes
it possible to:

* **unit-test the boundary algebra in isolation**, without running a
  full operator matvec on a full mesh;
* **swap closures** for ablation experiments
  (:class:`DDExtrapolation` for Phase A; the trivial
  :class:`CellCenter` strategy reproduces the pre-fix Defect-1
  behaviour for back-bisection);
* **localise the boundary algebra** so the matvec orchestration in
  :mod:`orpheus.sn.operator` stays thin and the closure-specific
  algebra lives next to the protocol that defines it.

Phase A scope — what this module does NOT touch
================================================

This module addresses Issue #168 Defects 1 + 2 (outer-face flux at
:math:`i = N - 1`, and the storage-conflation that corrupts the
interior-face stencil at :math:`i = N - 2` for inward ordinates).
Defect 3 (the sphere-pole stencil at :math:`i = 0`) is **preserved
as-is** in this phase, pending Phase B literature consultation
(Lewis & Miller §4.5; Bailey 2009; Carlson starting-direction;
Lathrop pole stencil).  The :meth:`DDExtrapolation.__call__` at
``cell_idx = 0`` returns the cell-centre value, which is correct at
the pole because the inner face area :math:`A_0` is zero by
construction in spherical geometry — the spatial flux at :math:`r = 0`
multiplies by zero.  The pole's Defect-3 mismatch is a separate
issue in the redistribution term, not in this face-flux closure.

Architectural mirror of CellUpdate
==================================

The Protocol shape mirrors
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` for consistency:

* ``@runtime_checkable Protocol`` so concrete strategies are matched
  by structural typing — third-party strategies / test mocks need
  not inherit.
* :class:`BoundaryFaceFluxBase` is a parallel concrete ABC layered
  on :class:`~orpheus.numerics.registry.RegistryMixin`; subclasses
  declare a ``key=`` class-creation kwarg and self-register, gaining
  name-keyed lookup via :meth:`BoundaryFaceFluxBase.create`.
* A class-level ``is_linear: bool`` trait advertises whether the
  closure is linear in the cell-centre array.  Both ``DDExtrapolation``
  and ``CellCenter`` are linear (affine combinations of cell-centre
  values).

References
==========

* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §5.3 (Diamond Difference closure).
* :doc:`/theory/discrete_ordinates` — "Boundary face-flux strategies
  for the curvilinear FD operator".
* Issue #168 — three independent boundary truncation defects;
  ``.claude/agent-memory/numerics-investigator/issue_168_three_defects.md``.
* Issue #168 design memo — ``.claude/plans/issue_168_design.md``.

See also
========

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the
  per-cell update strategy this module's contract mirrors.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar, Protocol, runtime_checkable

import numpy as np

from orpheus.numerics.registry import RegistryMixin


# ═══════════════════════════════════════════════════════════════════════
# BoundaryFaceFlux Protocol
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class BoundaryFaceFlux(Protocol):
    r"""Strategy contract for closing the radial face flux at a boundary cell.

    Concrete strategies (e.g. :class:`DDExtrapolation`) take the
    radial cell-centre array along the radial axis for one ordinate
    and one group, plus an integer ``cell_idx`` selecting either the
    inner pole (``cell_idx = 0``) or the outer face cell (``cell_idx
    = nx - 1``), and return the single face-flux value at the
    geometry boundary.

    The reduced ``(psi_cells, cell_idx)`` signature was chosen over
    the design-memo's five-argument ``(psi_cells, cell_idx, side,
    ord_idx, bc)`` because (a) ``side`` is derivable from
    ``cell_idx`` against a known ``nx``; (b) ``ord_idx`` was purely
    diagnostic and never read by the closure math; (c) ``bc`` is
    only consulted for the inner pole, which Phase A does not touch.
    The minimal signature gives the strategy exactly the data the
    DD-extrapolation algebra needs and nothing more.

    Attributes
    ----------
    is_linear : bool
        Whether the closure is linear (affine in the cell-centre
        array). DD extrapolation and the legacy cell-centre
        reproducer are both linear; future curvature-corrected
        closures could be non-linear (e.g. positivity-preserving
        clipped extrapolation).  Read-only class attribute.

    Notes
    -----
    The strategy is **stateless** — the same instance can be reused
    across every boundary call of a matvec without per-call
    allocation.  No ``LinearOperator`` capability tagging applies
    because a :class:`BoundaryFaceFlux` is a per-call scalar callable,
    not a vector-space operator with domain / range / inverse.
    """

    is_linear: bool

    def __call__(
        self,
        psi_cells: np.ndarray,
        cell_idx: int,
    ) -> np.ndarray:
        r"""Compute the radial face-flux value at this boundary cell.

        Parameters
        ----------
        psi_cells : np.ndarray
            Radial cell-centre angular flux for **one ordinate**, all
            groups.  Shape ``(ng, nx)`` — :math:`\psi_{n,i,g}` indexed
            by ``[g, i]``.
        cell_idx : int
            Which boundary cell.  ``0`` for the inner-pole (sphere
            r=0 / cylinder centerline) face; ``nx - 1`` for the outer
            face (geometry boundary).  Cells in between are not
            boundary cells — interior face fluxes are computed by the
            symmetric arithmetic average inside the matvec, not
            through this strategy.

        Returns
        -------
        np.ndarray
            Shape ``(ng,)`` — the per-group face flux at the
            specified boundary face.
        """
        ...


# ═══════════════════════════════════════════════════════════════════════
# BoundaryFaceFluxBase — concrete ABC with self-registration
# ═══════════════════════════════════════════════════════════════════════


class BoundaryFaceFluxBase(RegistryMixin, ABC):
    r"""Concrete abstract base for self-registering boundary-face-flux strategies.

    Subclasses inherit this ABC and pass ``key="..."`` in the class
    statement to self-register; the registry is consulted via
    :meth:`BoundaryFaceFluxBase.create("dd_extrapolation")` (or any
    other registered key).

    Subclasses MUST declare:

    * ``is_linear: ClassVar[bool]`` — whether the closure is linear
      in ``psi_cells``.
    * :meth:`__call__` — the per-boundary-cell algorithm.

    Notes
    -----
    Adding a new strategy is a one-line edit::

        class CurvatureCorrected(BoundaryFaceFluxBase, key="curvature_corrected"):
            is_linear: ClassVar[bool] = False

            def __call__(self, psi_cells, cell_idx):
                ...

    No registry insert; ``BoundaryFaceFluxBase.create(
    "curvature_corrected")`` is immediately callable.
    """

    registry: ClassVar[dict[str, type["BoundaryFaceFluxBase"]]] = {}

    is_linear: ClassVar[bool]

    @classmethod
    def _registry_base(cls) -> type:
        return BoundaryFaceFluxBase

    @abstractmethod
    def __call__(
        self,
        psi_cells: np.ndarray,
        cell_idx: int,
    ) -> np.ndarray:
        ...


# ═══════════════════════════════════════════════════════════════════════
# DDExtrapolation — Phase A default closure
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class DDExtrapolation(
    BoundaryFaceFluxBase, key="dd_extrapolation",
):
    r"""One-sided second-order Diamond-Difference boundary extrapolation.

    Phase A default for the curvilinear FD operator at the **outer
    face** (``cell_idx = nx - 1``):

    .. math::

       \psi^{\rm face}_{N-1/2} \;=\; \tfrac{3}{2}\,\psi_{N-1}
           \;-\; \tfrac{1}{2}\,\psi_{N-2}.

    This is the DD diamond relation
    :math:`\psi^{\rm face}_{i+1/2} = 2\,\bar\psi_i -
    \psi^{\rm face}_{i-1/2}` solved for the outer face when the
    interior face value :math:`\tfrac{1}{2}(\psi_{N-2} + \psi_{N-1})`
    is substituted for :math:`\psi^{\rm face}_{i-1/2}` — equivalent
    to a one-sided second-order finite-difference extrapolation
    through two cell centres. Truncation residual at the outer face
    is :math:`\mathcal{O}(h^2)`, matching the interior arithmetic
    average.

    At the **inner pole** (``cell_idx = 0``), Phase A preserves the
    historical Bailey treatment: return :math:`\psi_0`, the
    cell-centre value.  This is correct at :math:`r = 0` because the
    inner face area :math:`A_0` vanishes in spherical / cylindrical
    geometry — the spatial flux at the pole is multiplied by zero in
    the matvec, so its precise value is irrelevant.  Defect 3 (the
    sphere-pole redistribution-term mismatch) lives in a different
    term and is the subject of Phase B; this strategy intentionally
    does **NOT** address it.

    Notes
    -----
    Frozen + slotted: instances are immutable and lightweight; a
    single :class:`DDExtrapolation` can be reused across every
    boundary call of a matvec without per-call allocation.
    """

    is_linear: ClassVar[bool] = True
    """DD extrapolation is an affine combination of cell-centre
    values (coefficients :math:`\\tfrac{3}{2}, -\\tfrac{1}{2}`); the
    output depends linearly on ``psi_cells``."""

    def __call__(
        self,
        psi_cells: np.ndarray,
        cell_idx: int,
    ) -> np.ndarray:
        r"""Compute the boundary face flux via DD extrapolation.

        Parameters
        ----------
        psi_cells : np.ndarray
            Cell-centre angular flux along the radial axis, shape
            ``(ng, nx)``.
        cell_idx : int
            ``0`` for the inner pole (returns ``psi_cells[:, 0]`` as
            the Phase-A pole closure); ``nx - 1`` for the outer face
            (returns the DD extrapolation
            :math:`1.5\,\psi_{N-1} - 0.5\,\psi_{N-2}`).

        Returns
        -------
        np.ndarray
            Shape ``(ng,)``.
        """
        if cell_idx == 0:
            # Phase A — preserve the Bailey pole closure.  Defect 3
            # (sphere-pole redistribution-term mismatch) is deferred
            # to Phase B with literature consultation; the value
            # returned here is multiplied by A_0 = 0 in the spherical
            # matvec, so the closure is consistent with the geometry
            # even though the redistribution term itself is still
            # under the ERR-026 PARTIAL-CLOSURE flag.
            return psi_cells[:, 0]

        # Outer face — one-sided second-order extrapolation through
        # the two outermost cell centres.  Equivalent to the DD
        # diamond relation
        # ``psi_face_out = 2 * psi_cell - psi_face_in``
        # with the interior face value
        # ``psi_face_in = 0.5 * (psi[nx-2] + psi[nx-1])``
        # substituted in.
        return 1.5 * psi_cells[:, cell_idx] - 0.5 * psi_cells[:, cell_idx - 1]

    def __repr__(self) -> str:
        return "DDExtrapolation()"


# ═══════════════════════════════════════════════════════════════════════
# CellCenter — legacy reproducer for ablation / back-bisection
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class CellCenter(
    BoundaryFaceFluxBase, key="cell_center",
):
    r"""Legacy cell-centre closure — reproduces the Defect-1 behaviour.

    Returns ``psi_cells[:, cell_idx]`` for both inner pole and outer
    face.  This is the **pre-Phase-A** behaviour of the curvilinear
    FD operator — used the cell-centre as the face-flux value at
    :math:`i = N - 1` for outgoing :math:`\mu > 0`, giving
    :math:`\mathcal{O}(h)` truncation at the outer face.

    Why this strategy exists
    ------------------------

    * **Ablation studies**: comparing
      :class:`DDExtrapolation` vs :class:`CellCenter` quantifies the
      Phase-A improvement on convergence-rate at the boundary cell
      — the diagnostic-script mirror of the
      ``diag_issue168_*`` truncation-residual scaling tables.
    * **Back-bisection**: if a future change to the matvec causes a
      regression, swapping in :class:`CellCenter` reproduces the
      historical pre-Phase-A operator output and pins which step in
      a multi-commit refactor introduced the drift.
    * **Documentation**: the strategy carries the bug it reproduces
      in its docstring.  Issue #168 Defect 1 is now a strategy
      pluggable into the matvec — anyone reading
      :mod:`orpheus.sn.spatial.boundary_face_flux` learns what the
      bug was, not just what the fix is.

    This is **not** a default; production callers MUST use
    :class:`DDExtrapolation`.
    """

    is_linear: ClassVar[bool] = True
    """Identity on the cell-centre value is trivially linear."""

    def __call__(
        self,
        psi_cells: np.ndarray,
        cell_idx: int,
    ) -> np.ndarray:
        r"""Return the cell-centre value as a face-flux approximation.

        Parameters
        ----------
        psi_cells : np.ndarray
            Cell-centre angular flux along the radial axis, shape
            ``(ng, nx)``.
        cell_idx : int
            ``0`` or ``nx - 1`` (boundary cell).

        Returns
        -------
        np.ndarray
            ``psi_cells[:, cell_idx]``, shape ``(ng,)``.
        """
        return psi_cells[:, cell_idx]

    def __repr__(self) -> str:
        return "CellCenter()"


__all__ = [
    "BoundaryFaceFlux",
    "BoundaryFaceFluxBase",
    "CellCenter",
    "DDExtrapolation",
]

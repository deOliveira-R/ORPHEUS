r"""Angular flux at SN domain boundaries (and curvilinear pole state).

Issue #197 PR-TYPED-2 — the typed wrapper that retires the
``psi_bc: dict`` parameter previously threaded through
:func:`~orpheus.sn.sweep.transport_sweep` and stored on
:class:`~orpheus.sn.solver.SNSolver` as ``self._psi_bc``.  The dict
was a stringly-typed bag (``"bc_1d_left_face"``, ``"bc_2d_x"``,
``"psi_pole"``, ``"phi_0_prev"``, ...): a typo like
``"bc_xmim"`` (instead of ``"bc_xmin"``) silently lazy-initialised
a fresh zero buffer at the wrong key, losing the persistent-BC
state without surfacing an error.  Per ``coding-elegance`` Pattern 4
(illegal states unrepresentable), the dict gives way to named
fields on a dataclass; typos become AttributeErrors at write time.

Shape contracts
===============

The per-face shapes match the sweep's working layout (Issue #196
PR-INDEX-1..5/7).

* **1-D slab** — two faces, ``xmin`` and ``xmax``, each ``(N, ng)``.
* **1-D spherical / cylindrical** — one outer radial face ``xmax``
  ``(N, ng)``.  Curvilinear pole iteration history lives on
  :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
  (R-1 Step 0, 2026-05-19) — see the module docstring.
* **2-D Cartesian** — two persistent buffers covering BOTH x-faces
  AND interior x-face cells (``xmin_xmax_buf`` shape ``(N, ng, nx+1, ny)``)
  and likewise for y (``ymin_ymax_buf`` shape ``(N, ng, nx, ny+1)``).
  The ``xmin`` / ``xmax`` / ``ymin`` / ``ymax`` face *slices* are
  exposed via the accessor properties below; the persistent storage
  is the full buffer because the 2-D wavefront sweep writes
  outgoing-face cells back into it for the next reflective-BC apply.

Mutability
==========

:class:`BoundaryFlux` is **mutable** by design: the sweep's persistent
BC contract is a write-through cache.  Outgoing-face cells are written
back into the same buffers between calls, so reflective BC partners
read the most-recent outgoing flux on each subsequent sweep.
Wrapping the buffers in a frozen dataclass would force every sweep to
allocate fresh buffers and discard them — a memory churn the
production hot path cannot afford.

R-1 Step 0 (2026-05-19) — the curvilinear pole iteration history that
formerly lived here as ``pole_psi`` / ``pole_phi_prev`` has been moved to
:class:`~orpheus.sn.spatial.pole_angular_closure._PoleHistoryState` on
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`.
That state was always a curvilinear-closure concept; its presence on
:class:`BoundaryFlux` was a concept leak between two semantically distinct
notions (boundary face flux vs pole iteration cache).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh


__all__ = ["BoundaryFlux"]


@dataclass
class BoundaryFlux:
    r"""Angular flux at SN boundaries plus curvilinear pole state.

    Parameters
    ----------
    mesh : SNMesh
        The SN phase-space carrier.

    Attributes
    ----------
    xmin_face, xmax_face : np.ndarray or None
        1-D face buffers, shape ``(N, ng)``.  Populated for 1-D meshes
        (slab, spherical, cylindrical).  ``None`` for 2-D meshes (the
        2-D x-faces live inside :attr:`xmin_xmax_buf`).
    xmin_xmax_buf, ymin_ymax_buf : np.ndarray or None
        2-D persistent buffers covering both the face AND the interior
        face-flux cache the wavefront sweep maintains.  Shapes
        ``(N, ng, nx+1, ny)`` and ``(N, ng, nx, ny+1)``.  ``None`` for
        1-D meshes.
    """

    mesh: "SNMesh"

    # 1-D face buffers (slab + curvilinear outer radial face).
    xmin_face: np.ndarray | None = None
    xmax_face: np.ndarray | None = None

    # 2-D persistent buffers (cover face + interior face-flux cache).
    xmin_xmax_buf: np.ndarray | None = None
    ymin_ymax_buf: np.ndarray | None = None

    # ── Construction helpers ──────────────────────────────────────────

    @classmethod
    def zeros(cls, mesh: "SNMesh") -> "BoundaryFlux":
        r"""Build an all-zeros :class:`BoundaryFlux` sized to ``mesh``.

        Allocates only the buffers the mesh's geometry actually uses:

        * 1-D Cartesian (slab): ``xmin_face`` + ``xmax_face`` only.
        * 1-D spherical / cylindrical: ``xmax_face`` only (the
          ``xmin_face`` slot is unused — the geometric pole is not a
          BC face, it's a regularity condition).
        * 2-D Cartesian: ``xmin_xmax_buf`` + ``ymin_ymax_buf``.

        Pole state stays ``None`` until the first sweep writes it.
        """
        N = mesh.quad.N
        ng = mesh.ng
        nx = mesh.nx
        ny = mesh.ny
        reduced = mesh.reduced
        bf = cls(mesh=mesh)
        if reduced is not None:
            # 1-D — slab or curvilinear.
            curv = getattr(mesh, "curvature", None)
            if curv in ("spherical", "cylindrical"):
                # Curvilinear: only the outer radial face exists.
                bf.xmax_face = np.zeros((N, ng))
            else:
                # 1-D slab: two faces.
                bf.xmin_face = np.zeros((N, ng))
                bf.xmax_face = np.zeros((N, ng))
        else:
            # 2-D Cartesian: persistent buffers covering face + interior.
            bf.xmin_xmax_buf = np.zeros((N, ng, nx + 1, ny))
            bf.ymin_ymax_buf = np.zeros((N, ng, nx, ny + 1))
        return bf

    # ── Accessors (face slices for the 2-D persistent buffers) ────────
    #
    # The 2-D persistent buffers carry interior cells too; ``xmin`` /
    # ``xmax`` / ``ymin`` / ``ymax`` here expose the FACE slice that
    # the BC operator consumes.  They return ``np.ndarray`` views, not
    # owning copies — assignment to these slices propagates to the
    # underlying buffer (the intended write-through path).
    #
    # For 1-D meshes, ``xmin`` / ``xmax`` route to :attr:`xmin_face` /
    # :attr:`xmax_face` directly; 1-D meshes have no ``ymin`` / ``ymax``.

    @property
    def xmin(self) -> np.ndarray | None:
        r""":math:`\psi` on the ``xmin`` face.

        1-D: the ``(N, ng)`` :attr:`xmin_face` buffer (slab only).
        2-D: a view of :attr:`xmin_xmax_buf` at face slot ``[:, :, 0, :]``
        with shape ``(N, ng, ny)``.
        """
        if self.xmin_xmax_buf is not None:
            return self.xmin_xmax_buf[:, :, 0, :]
        return self.xmin_face

    @property
    def xmax(self) -> np.ndarray | None:
        r""":math:`\psi` on the ``xmax`` face.

        1-D: the ``(N, ng)`` :attr:`xmax_face` buffer (slab + curvilinear).
        2-D: a view of :attr:`xmin_xmax_buf` at face slot ``[:, :, nx, :]``
        with shape ``(N, ng, ny)``.
        """
        if self.xmin_xmax_buf is not None:
            nx = self.mesh.nx
            return self.xmin_xmax_buf[:, :, nx, :]
        return self.xmax_face

    @property
    def ymin(self) -> np.ndarray | None:
        r""":math:`\psi` on the ``ymin`` face (2-D only).

        A view of :attr:`ymin_ymax_buf` at face slot ``[:, :, :, 0]``
        with shape ``(N, ng, nx)``.  ``None`` for 1-D meshes.
        """
        if self.ymin_ymax_buf is not None:
            return self.ymin_ymax_buf[:, :, :, 0]
        return None

    @property
    def ymax(self) -> np.ndarray | None:
        r""":math:`\psi` on the ``ymax`` face (2-D only).

        A view of :attr:`ymin_ymax_buf` at face slot ``[:, :, :, ny]``
        with shape ``(N, ng, nx)``.  ``None`` for 1-D meshes.
        """
        if self.ymin_ymax_buf is not None:
            ny = self.mesh.ny
            return self.ymin_ymax_buf[:, :, :, ny]
        return None

    # ── R-1 Step 1 — algebra-element arithmetic ──────────────────────
    #
    # :class:`BoundaryFlux` participates in the operator algebra as the
    # face-flux companion to :class:`~orpheus.sn.angular_flux.AngularFlux`.
    # Arithmetic on :class:`AngularFlux` propagates to its ``.boundary``
    # field; that propagation delegates to the dunders below.
    #
    # Semantics: elementwise on every non-None face buffer.  The result
    # is a fresh :class:`BoundaryFlux` (immutable arithmetic).  The
    # sweep's mutable write-through pattern still works because it
    # operates on whichever instance it's handed — the arithmetic just
    # produces a new instance.
    #
    # Mesh-identity is the validity invariant: two operands must share
    # the SAME :class:`SNMesh` object (``is`` comparison) to combine.
    # This matches the existing :class:`AngularFlux._validate_partner`
    # pattern and ensures shape compatibility by construction.

    def _validate_partner(self, other: "BoundaryFlux") -> None:
        """Reject cross-mesh arithmetic.

        Two :class:`BoundaryFlux` instances may be combined only when
        they share the same :class:`SNMesh` reference (identity ``is``).
        Different meshes carry different geometry / quadrature / size,
        so combining them is meaningless.  Matches
        :meth:`AngularFlux._validate_partner`.
        """
        if not isinstance(other, BoundaryFlux):
            raise TypeError(
                f"BoundaryFlux arithmetic operand must be BoundaryFlux; "
                f"got {type(other).__name__}"
            )
        if self.mesh is not other.mesh:
            raise ValueError(
                "BoundaryFlux operands must share the same SNMesh "
                "reference (mesh-identity invariant); operating on "
                "BoundaryFlux instances bound to different meshes is "
                "undefined."
            )

    def _binary_op(
        self,
        other: "BoundaryFlux",
        op,
    ) -> "BoundaryFlux":
        """Element-wise binary op across every non-None face buffer."""
        self._validate_partner(other)
        return BoundaryFlux(
            mesh=self.mesh,
            xmin_face=(
                op(self.xmin_face, other.xmin_face)
                if self.xmin_face is not None and other.xmin_face is not None
                else None
            ),
            xmax_face=(
                op(self.xmax_face, other.xmax_face)
                if self.xmax_face is not None and other.xmax_face is not None
                else None
            ),
            xmin_xmax_buf=(
                op(self.xmin_xmax_buf, other.xmin_xmax_buf)
                if self.xmin_xmax_buf is not None and other.xmin_xmax_buf is not None
                else None
            ),
            ymin_ymax_buf=(
                op(self.ymin_ymax_buf, other.ymin_ymax_buf)
                if self.ymin_ymax_buf is not None and other.ymin_ymax_buf is not None
                else None
            ),
        )

    def _scalar_op(self, scalar: float, op) -> "BoundaryFlux":
        """Element-wise scalar op across every non-None face buffer."""
        return BoundaryFlux(
            mesh=self.mesh,
            xmin_face=op(self.xmin_face, scalar) if self.xmin_face is not None else None,
            xmax_face=op(self.xmax_face, scalar) if self.xmax_face is not None else None,
            xmin_xmax_buf=(
                op(self.xmin_xmax_buf, scalar)
                if self.xmin_xmax_buf is not None else None
            ),
            ymin_ymax_buf=(
                op(self.ymin_ymax_buf, scalar)
                if self.ymin_ymax_buf is not None else None
            ),
        )

    def __add__(self, other: "BoundaryFlux") -> "BoundaryFlux":
        return self._binary_op(other, lambda a, b: a + b)

    def __sub__(self, other: "BoundaryFlux") -> "BoundaryFlux":
        return self._binary_op(other, lambda a, b: a - b)

    def __mul__(self, scalar: float) -> "BoundaryFlux":
        return self._scalar_op(scalar, lambda a, c: a * c)

    def __rmul__(self, scalar: float) -> "BoundaryFlux":
        return self._scalar_op(scalar, lambda a, c: c * a)

    def __truediv__(self, scalar: float) -> "BoundaryFlux":
        return self._scalar_op(scalar, lambda a, c: a / c)

    def __neg__(self) -> "BoundaryFlux":
        return self._scalar_op(0.0, lambda a, _c: -a)

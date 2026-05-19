r"""Angular flux field on an SN phase space.

Issue #197 PR-TYPED-2 — the typed wrapper for
:math:`\psi(\vec r, \hat\Omega_n, g)` sampled on the SN phase space
(quadrature × spatial grid × energy groups).

R-1 Step 2 (2026-05-19) — :class:`AngularFlux` now **carries its own
boundary state** via a typed ``boundary: BoundaryFlux`` field.  The
:class:`BoundaryFlux` is auto-allocated as zeros at construction when
not supplied (``__post_init__`` invokes
:meth:`BoundaryFlux.zeros` on the mesh).  Arithmetic dunders propagate
to both ``.values`` and ``.boundary`` element-wise — the operator
algebra reads as math (``L.apply(psi)`` consumes the typed pair
without a second argument).

Frozen dataclass (``coding-elegance`` Pattern 4 — illegal states
unrepresentable): every :class:`AngularFlux` instance has a non-None
``boundary`` by the time anything reads it; mesh-shape mismatches
surface at construction time.

Reads as the math via dunder arithmetic
(``coding-elegance`` Pattern 1 — match the algebra of the domain) and
a typed angular reduction (``coding-elegance`` Pattern 3 — named
intermediates):

* ``psi_a + psi_b`` returns a new :class:`AngularFlux` with both
  cell-centre AND boundary arithmetic propagated.
* ``alpha * psi`` returns a new :class:`AngularFlux`.
* ``psi.integrate_angular()`` reduces along the ordinate axis,
  returning a :class:`ScalarFlux` (:math:`\phi = \sum_n w_n \psi_n`).
* ``psi.at_ordinate(n)`` returns the per-ordinate ``(ng, nx, ny)``
  slice view.

Flat-vector adapters for scipy.gmres integration:

* :meth:`AngularFlux.from_flat_with_traces` decodes a B1''-aware
  packed 1-D vector (cell + outer-face + inner-face blocks) into an
  :class:`AngularFlux` with face state in ``boundary``.
* :meth:`AngularFlux.to_flat_with_traces` is the inverse — concatenates
  cell values + boundary face values into the flat layout that
  scipy.gmres consumes.

These two methods are the ONLY ravel/reshape sites in the architecture
— operator-algebra consumers see typed :class:`AngularFlux` end-to-end.

Units: :math:`[1/(\rm cm^2\,s\,sr\,eV)]` per energy group bin.  The
per-bin energy density is absorbed into the cross-section convention.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .boundary_flux import BoundaryFlux
    from .geometry import SNMesh
    from .scalar_flux import ScalarFlux


__all__ = ["AngularFlux"]


@dataclass(frozen=True)
class AngularFlux:
    r"""Angular flux field :math:`\psi(\vec r, \hat\Omega_n, g)`.

    Parameters
    ----------
    values : np.ndarray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    mesh : SNMesh
        The SN phase-space carrier — validates shape and supplies the
        canonical ``(N, ng, nx, ny)`` quadruple.
    boundary : BoundaryFlux or None, optional
        Boundary face flux companion.  When ``None`` (default), a
        fresh :class:`BoundaryFlux.zeros(mesh)` is allocated by
        ``__post_init__``.  Pattern 4: every :class:`AngularFlux`
        instance HAS a boundary by the time anything reads it.
    """

    values: np.ndarray
    mesh: "SNMesh"
    boundary: "BoundaryFlux | None" = None

    def __post_init__(self) -> None:
        # Shape invariant.
        N = self.mesh.quad.N
        expected = (N, self.mesh.ng, self.mesh.nx, self.mesh.ny)
        if self.values.shape != expected:
            raise ValueError(
                f"AngularFlux expects shape (N, ng, nx, ny) = {expected}; "
                f"got {self.values.shape}"
            )
        # R-1 Step 2: auto-allocate boundary when not supplied.
        # ``object.__setattr__`` bypasses frozen-dataclass restriction
        # (same pattern Mesh1D uses for its eager attributes).
        if self.boundary is None:
            from .boundary_flux import BoundaryFlux
            object.__setattr__(
                self, "boundary", BoundaryFlux.zeros(self.mesh),
            )

    # ── Reduction → ScalarFlux ────────────────────────────────────────

    def integrate_angular(self) -> "ScalarFlux":
        r"""Reduce along the ordinate axis: :math:`\phi = \sum_n w_n \psi_n`.

        Returns a :class:`~orpheus.sn.scalar_flux.ScalarFlux` of shape
        ``(ng, nx, ny)``.  The per-ordinate quadrature weights live on
        the mesh's :class:`AngularQuadrature`; the contraction is the
        canonical angular average that produces the scalar flux that
        the within-group scattering source consumes.
        """
        from .scalar_flux import ScalarFlux
        w = self.mesh.quad.weights
        return ScalarFlux(
            np.einsum("n,ngxy->gxy", w, self.values), self.mesh,
        )

    # ── Dunder algebra (within-type) ──────────────────────────────────
    #
    # R-1 Step 2: arithmetic propagates to ``boundary`` by delegating to
    # :class:`BoundaryFlux`'s dunders (R-1 Step 1).  The pair
    # (values, boundary) is the algebraic state; the propagation rule
    # is one line.

    def __add__(self, other: "AngularFlux") -> "AngularFlux":
        self._validate_partner(other)
        return AngularFlux(
            values=self.values + other.values,
            mesh=self.mesh,
            boundary=self.boundary + other.boundary,
        )

    def __sub__(self, other: "AngularFlux") -> "AngularFlux":
        self._validate_partner(other)
        return AngularFlux(
            values=self.values - other.values,
            mesh=self.mesh,
            boundary=self.boundary - other.boundary,
        )

    def __mul__(self, scalar: float) -> "AngularFlux":
        c = float(scalar)
        return AngularFlux(
            values=self.values * c,
            mesh=self.mesh,
            boundary=self.boundary * c,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "AngularFlux":
        c = float(scalar)
        return AngularFlux(
            values=self.values / c,
            mesh=self.mesh,
            boundary=self.boundary / c,
        )

    def __neg__(self) -> "AngularFlux":
        return AngularFlux(
            values=-self.values,
            mesh=self.mesh,
            boundary=-self.boundary,
        )

    def _validate_partner(self, other: "AngularFlux") -> None:
        if not isinstance(other, AngularFlux):
            raise TypeError(
                "AngularFlux arithmetic requires an AngularFlux partner; "
                f"got {type(other).__name__}"
            )
        if other.mesh is not self.mesh:
            raise ValueError(
                "AngularFlux arithmetic across distinct SNMesh instances "
                "is forbidden — the field is mesh-bound."
            )

    # ── Flat-vector adapters (scipy boundary) ─────────────────────────
    #
    # The B1''-aware methods are the ONLY ravel/reshape sites in the
    # entire architecture.  Operator-algebra consumers see typed
    # :class:`AngularFlux` end-to-end; only the scipy.gmres adapter (or
    # similar dense-linear-algebra integration) flattens at the boundary.
    #
    # Two flavours for the two extant equation-map layouts:
    #
    # * ``_with_traces`` — B1'' face-aware layout used by
    #   :class:`~orpheus.sn.operator.StreamingOperator` /
    #   :class:`~orpheus.sn.operator.CollisionOperator` on 1-D meshes.
    #   ``n_unknowns = (N·nx + n_face_outer + n_face_inner) · ng``.
    # * ``_legacy`` — pre-B1'' cell-only layout used by
    #   :class:`~orpheus.sn.operator.SNStreamingOperator` and the 2-D
    #   Cartesian wavefront path.  ``n_unknowns = n_eq · ng`` where
    #   ``n_eq`` may be COMPRESSED on curvilinear (drops inward-at-outer
    #   slots that the analytical extension fills).  Each layout has its
    #   own factory name (Pattern 4 — illegal states unrepresentable).

    @classmethod
    def from_flat_with_traces(
        cls, flat: np.ndarray, mesh: "SNMesh",
    ) -> "AngularFlux":
        r"""Decode a B1''-aware packed 1-D vector into an :class:`AngularFlux`.

        The packed layout contains three concatenated blocks:

        1. Cell-centre block: ``(N·nx · ng,)`` — every ``(ordinate, cell)``
           pair is an unknown.
        2. Outer-face block: ``(n_face_outer · ng,)`` — one slot per
           outward ordinate at the outer boundary (``r=R`` / ``x=L``).
        3. Inner-face block (slab only): ``(n_face_inner · ng,)`` — one
           slot per inward ordinate at the inner boundary (``x=0``).

        Decoded into a fresh :class:`AngularFlux` whose ``.boundary``
        carries the face state at ``xmax_face`` (outer) and (slab only)
        ``xmin_face`` (inner).

        Parameters
        ----------
        flat
            Packed 1-D vector matching the B1'' layout for ``mesh``.
        mesh
            The :class:`SNMesh` whose geometry determines the eq_map
            (and therefore the expected ``flat.size``).

        Returns
        -------
        AngularFlux
            Typed angular flux with B1'' face state in ``boundary``.
        """
        from .boundary_flux import BoundaryFlux
        from .operator import (
            build_equation_map_with_traces,
            solution_to_angular_flux_with_traces,
        )
        nx, ng = mesh.nx, mesh.ng
        N = mesh.quad.N
        has_inner_bc = getattr(mesh, "curvature", None) is None
        eq_map = build_equation_map_with_traces(
            nx, mesh.quad, ng, has_inner_bc=has_inner_bc,
        )
        if flat.size != eq_map.n_unknowns:
            raise ValueError(
                f"AngularFlux.from_flat_with_traces: flat.size = {flat.size} "
                f"does not match eq_map.n_unknowns = {eq_map.n_unknowns} "
                f"for mesh (N={N}, nx={nx}, ng={ng}, "
                f"has_inner_bc={has_inner_bc})."
            )
        psi_cell, psi_face_outer, psi_face_inner = (
            solution_to_angular_flux_with_traces(
                flat, eq_map, nx, ng, N=N,
            )
        )
        # Scatter face values into BoundaryFlux face buffers.
        boundary = BoundaryFlux.zeros(mesh)
        if eq_map.n_face_outer > 0:
            boundary.xmax_face[eq_map.face_outer_ordinate, :] = psi_face_outer
        if eq_map.n_face_inner > 0 and psi_face_inner is not None:
            boundary.xmin_face[eq_map.face_inner_ordinate, :] = psi_face_inner
        return cls(values=psi_cell, mesh=mesh, boundary=boundary)

    def to_flat_with_traces(self) -> np.ndarray:
        r"""Encode this :class:`AngularFlux` into a B1''-aware packed 1-D vector.

        Inverse of :meth:`from_flat_with_traces`: concatenates the cell
        block (gathered at the eq_map's ``(ordinate, ix)`` slots) with
        the outer-face and (slab only) inner-face blocks (gathered at
        the per-face-ordinate slots).

        Round-trip identity (pinned in
        :file:`tests/sn/test_angular_flux_with_boundary.py`):
        ``AngularFlux.from_flat_with_traces(flat, mesh).to_flat_with_traces() == flat``.
        """
        from .operator import (
            build_equation_map_with_traces,
            pack_with_traces,
        )
        nx, ng = self.mesh.nx, self.mesh.ng
        has_inner_bc = getattr(self.mesh, "curvature", None) is None
        eq_map = build_equation_map_with_traces(
            nx, self.mesh.quad, ng, has_inner_bc=has_inner_bc,
        )
        # Gather face values from the boundary buffers at the typed
        # per-ordinate slots.  Face buffer indexing relies on the
        # boundary already having ``xmax_face`` / (slab) ``xmin_face``
        # allocated — guaranteed by ``__post_init__``.
        if eq_map.n_face_outer > 0:
            face_outer = self.boundary.xmax_face[
                eq_map.face_outer_ordinate, :
            ]
        else:
            face_outer = np.zeros((0, ng))
        if eq_map.n_face_inner > 0:
            face_inner = self.boundary.xmin_face[
                eq_map.face_inner_ordinate, :
            ]
        else:
            face_inner = None
        return pack_with_traces(
            self.values, face_outer, face_inner, eq_map,
        )

    # ── Selectors ─────────────────────────────────────────────────────

    def at_ordinate(self, n: int) -> np.ndarray:
        """Return the per-ordinate slice ``values[n]``, shape ``(ng, nx, ny)``."""
        return self.values[n]

    # ── Diagnostics ───────────────────────────────────────────────────

    def linf(self) -> float:
        r"""Return :math:`\lVert\psi\rVert_\infty` over all entries."""
        return float(np.abs(self.values).max())

    def copy(self) -> "AngularFlux":
        """Return a deep copy carrying owned ``values`` AND ``boundary``.

        Both the cell-centre values and the boundary face state are
        cloned independently — the resulting :class:`AngularFlux` shares
        no memory with the original.  Useful for SI outer loops that
        carry the previous iterate as ``initial_guess`` without risk of
        the current sweep mutating it.
        """
        # ``boundary.copy()`` doesn't exist on BoundaryFlux (it's a
        # @dataclass; explicit ndarray copy).  Build a fresh BoundaryFlux
        # by cloning each non-None face buffer.
        from .boundary_flux import BoundaryFlux
        b = self.boundary
        new_b = BoundaryFlux(
            mesh=self.mesh,
            xmin_face=b.xmin_face.copy() if b.xmin_face is not None else None,
            xmax_face=b.xmax_face.copy() if b.xmax_face is not None else None,
            xmin_xmax_buf=(
                b.xmin_xmax_buf.copy() if b.xmin_xmax_buf is not None else None
            ),
            ymin_ymax_buf=(
                b.ymin_ymax_buf.copy() if b.ymin_ymax_buf is not None else None
            ),
        )
        return AngularFlux(self.values.copy(), self.mesh, boundary=new_b)

    # ── Metadata read-throughs ────────────────────────────────────────

    @property
    def N(self) -> int:  # noqa: N802 — matches AngularQuadrature.N
        """Number of angular ordinates."""
        return self.mesh.quad.N

    @property
    def ng(self) -> int:
        """Energy group count."""
        return self.mesh.ng

    @property
    def nx(self) -> int:
        """Spatial extent in x."""
        return self.mesh.nx

    @property
    def ny(self) -> int:
        """Spatial extent in y."""
        return self.mesh.ny

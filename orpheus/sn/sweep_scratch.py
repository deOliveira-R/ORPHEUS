r"""Sweep-private interior wavefront cache for the 2-D Cartesian sweep.

L3 (SN-method-specific). The 2-D wavefront sweep
(:func:`~orpheus.sn.sweep._sweep_2d_wavefront`) propagates angular
fluxes through cells along anti-diagonal levels. As it sweeps each
octant, it writes the per-cell face fluxes into a working buffer of
shape ``(N, ng, nx+1, ny)`` (for x-faces) and ``(N, ng, nx, ny+1)``
(for y-faces). The buffer slots ``[..., 0, :]`` / ``[..., nx, :]`` are
the boundary face slots (consumed by the reflective BC partner on the
next iteration); the slots ``[..., 1:nx, :]`` and ``[..., :, 1:ny]``
are the interior wavefront cache (sweep-private).

Motivation (Depth B step D-G)
==============================

Pre-D-G these working buffers lived on
:class:`~orpheus.sn.boundary_flux.BoundaryFlux` as
``xmin_xmax_buf`` / ``ymin_ymax_buf``. This conflated TWO concepts:

* **Boundary face trace** (sweep -> BoundaryFlux -> reflective BC) —
  load-bearing inter-iteration state.
* **Interior wavefront cache** (sweep -> sweep) — sweep-private
  scratch that does not need to persist across iterations as
  authoritative state (each sweep recomputes every interior cell from
  upstream propagation).

Conflating them forced BoundaryFlux to be mutable (because the sweep
needed write-through into the interior buffer), even though the
boundary-face concept is naturally immutable per-iteration. Pre-D-G
this constraint blocked making BoundaryFlux a pure
:class:`~orpheus.numerics.field.Field`.

Post-D-G :class:`SweepScratch` owns the FULL working buffers; the
sweep allocates a fresh :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
(face-only, small) on each iteration carrying the new face-slot values
that the BC partner reads next. SweepScratch persists across sweep
calls (owned by the sweep operator, not allocated per-call) — the
memory-churn constraint that pre-D-G mutability documented at
``boundary_flux.py`` lines 34-43 is preserved.

References
==========

* Depth B plan §3.4 (storage split).
* Plan §6 step D-G (BoundaryFlux as pure Field + SweepScratch carve).
* Pre-D-G :class:`~orpheus.sn.boundary_flux.BoundaryFlux` docstring
  lines 34-43 (mutability rationale).
* ``coding-elegance`` Pattern 3 (named intermediates with domain
  semantics) — boundary trace and interior cache are separately
  named after the carve.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .geometry import SNMesh


__all__ = ["SweepScratch"]


@dataclass
class SweepScratch:
    r"""Sweep-private interior wavefront cache for the 2-D Cartesian sweep.

    Attributes
    ----------
    psi_x_buf : np.ndarray or None
        Working buffer for x-direction face fluxes during a 2-D sweep,
        shape ``(N, ng, nx+1, ny)``. Slots ``[..., 0, :]`` and
        ``[..., nx, :]`` are the xmin / xmax boundary face slots
        (read-in from BoundaryFlux; written-out to a fresh BoundaryFlux
        at sweep end). Slots ``[..., 1:nx, :]`` are sweep-private
        interior wavefront cache (persist across calls in this scratch).
        ``None`` for 1-D meshes (the 1-D sweep does not use a 2-D
        wavefront buffer).
    psi_y_buf : np.ndarray or None
        Working buffer for y-direction face fluxes during a 2-D sweep,
        shape ``(N, ng, nx, ny+1)``. Slots ``[..., :, 0]`` and
        ``[..., :, ny]`` are ymin / ymax boundary face slots; slots
        ``[..., :, 1:ny]`` are sweep-private interior wavefront cache.
        ``None`` for 1-D meshes.

    Notes
    -----
    :class:`SweepScratch` is **mutable** by design and intended to be
    owned by the sweep operator across multiple sweep invocations
    (preserves the pre-D-G "no per-sweep large-buffer allocation"
    invariant — see module docstring). The pre-D-G mutability
    rationale at ``boundary_flux.py`` lines 34-43 applies here verbatim,
    but now constrained to the interior scratch only; BoundaryFlux is
    free to become a pure immutable :class:`~orpheus.numerics.field.Field`.

    Allocation policy: the :meth:`for_sn_mesh` classmethod allocates
    eagerly for 2-D meshes and returns empty (``None`` fields) for 1-D
    meshes. The sweep may lazy-init these fields if it receives an
    empty scratch for a 2-D mesh (matches the pre-D-G lazy-init pattern
    in ``sweep.py`` lines 769-772).
    """

    psi_x_buf: NDArray | None = None
    psi_y_buf: NDArray | None = None

    @classmethod
    def for_sn_mesh(cls, mesh: "SNMesh") -> "SweepScratch":
        r"""Allocate scratch sized to ``mesh``.

        Returns an empty :class:`SweepScratch` (``psi_x_buf`` and
        ``psi_y_buf`` both ``None``) for 1-D meshes (slab + curvilinear
        sphere / cylinder). Returns an allocated :class:`SweepScratch`
        with zero-initialised ``(N, ng, nx+1, ny)`` and
        ``(N, ng, nx, ny+1)`` buffers for 2-D Cartesian meshes.

        Parameters
        ----------
        mesh : SNMesh
            The SN phase-space carrier. Discriminates 1-D vs 2-D via
            ``mesh.reduced is None`` (2-D) or not None (1-D).

        Returns
        -------
        SweepScratch
            Scratch sized to ``mesh``.
        """
        if mesh.reduced is not None:
            # 1-D — no 2-D wavefront buffer needed.
            return cls()
        N = mesh.quad.N
        ng = mesh.ng
        nx = mesh.nx
        ny = mesh.ny
        return cls(
            psi_x_buf=np.zeros((N, ng, nx + 1, ny)),
            psi_y_buf=np.zeros((N, ng, nx, ny + 1)),
        )

    def ensure_2d_buffers(self, mesh: "SNMesh") -> None:
        r"""Lazy-init the 2-D buffers if not yet allocated.

        Matches the pre-D-G ``sweep.py`` lines 769-772 pattern (lazy-init
        on first sweep call). Idempotent — if buffers exist, no-op. If
        they don't and the mesh is 1-D, raises (1-D never uses these
        buffers).

        Parameters
        ----------
        mesh : SNMesh
            The SN phase-space carrier; used for shape determination.

        Raises
        ------
        ValueError
            If ``mesh.reduced is not None`` (1-D mesh) but the caller
            asked to allocate 2-D buffers.
        """
        if mesh.reduced is not None:
            raise ValueError(
                "SweepScratch.ensure_2d_buffers: requested 2-D buffer "
                "allocation for a 1-D mesh (mesh.reduced is not None)"
            )
        N = mesh.quad.N
        ng = mesh.ng
        nx = mesh.nx
        ny = mesh.ny
        if self.psi_x_buf is None:
            self.psi_x_buf = np.zeros((N, ng, nx + 1, ny))
        if self.psi_y_buf is None:
            self.psi_y_buf = np.zeros((N, ng, nx, ny + 1))

r"""Per-octant causal cell DAG for the 2-D Cartesian wavefront sweep.

This module ships the **§15A.2 "upwind trace complex / causal transport
DAG / direction sweep ordering"** primitive (Grand Report v3 lines
2137-2171) for 2-D Cartesian SN, lifting the per-call ``_diag_cache``
precompute that previously lived inside
:func:`orpheus.sn.sweep._sweep_2d_wavefront` (lines 766-785) to
mesh-time work, and replacing the per-ordinate Python loop over
ordinates with a per-octant batched ``apply``.

The graph is a **derived object** — it depends only on cell topology
(``Mesh2D``) and the octant sign convention. It does NOT depend on
fluxes, sources, or iteration state, so it is built ONCE per
``(SNMesh, octant)`` pair at mesh construction (Wave 2 / C2.4) and
reused across every source iteration / Krylov matvec / outer iteration.

Architectural framing (Cardinal Rule 2)
=======================================

Per the project memory note ``project_moc_structure.md`` and
:doc:`../sn/spatial/cell_update`, this DAG abstraction is **SN-specific
by design**. MoC will define its own analog (per-ray traversal — fiber
bundles + solution sheaves, NOT a topological sort over a cell graph).
There is no shared ``SweepGraph`` Protocol because there is no shared
mathematical structure.

Tensorial framing (Wave 2 plan, §15A.2)
=======================================

Per octant, the streaming inverse :math:`L^{-1}_{\text{oct}}` acts on
the ``(octant_ordinates × cells × groups)`` tensor::

    L⁻¹_oct : (N_oct, nx, ny, ng) -> (N_oct, nx, ny, ng)

The anti-diagonal schedule vectorises ALL THREE axes inside one
``update_batch`` call — no Python loop over ordinates, no Python loop
over cells along a diagonal, no Python loop over groups. The outer
sweep iterates octants (4 in 2-D — structural) and topological levels
per octant (structural — sweep-DAG data dependency). The ordinate
axis is internal to every ``apply`` / ``einsum`` invocation.

References
==========

* Grand Report v3 §15A.2 lines 2137-2171 — "upwind trace complex /
  causal transport DAG / direction sweep ordering" primitive
  description with the ``assert_*`` invariant set this module's
  tests pin (``assert_upwind_orientation``,
  ``assert_topologically_sorted``, ``assert_face_pairing_consistent``,
  ``assert_boundary_trace_classification``,
  ``assert_cycles_are_declared``).
* Wave 2 plan ``.claude/plans/transient-giggling-cake.md`` C2.3 —
  this module's design + dispatch boundary with ``CellUpdate``.
* :class:`~orpheus.sn.spatial.cell_update.SweepCellSlice` — the
  per-level packet ``apply`` builds.
* :meth:`~orpheus.sn.spatial.cell_update.CellUpdateBase.update_batch`
  — the closure-pluggable consumer (Wave 2 / C2.2).

See also
========

* Future direction (out of scope): assemble per-octant lower-triangular
  ``L_oct`` as ``scipy.sparse.csr_matrix`` and replace the
  forward-substitution ``apply`` with
  ``scipy.sparse.linalg.spsolve_triangular``. Trades memory for a
  compiled inner kernel; worth it for very large 3-D problems. File a
  follow-up when the 3-D consumer arrives.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from orpheus.sn.spatial.cell_update import CellUpdateBase, SweepCellSlice


# ═══════════════════════════════════════════════════════════════════════
# OctantLabel — 2-D in-plane octant signature
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class OctantLabel:
    r"""Octant signature for the 2-D Cartesian sweep.

    Carries the in-plane direction signs ``(sign_x, sign_y) ∈
    {-1, 0, +1}²``. The ``sign_z`` component of an ``S²`` ordinate
    label is dropped — the 2-D Cartesian sweep is invariant under
    the out-of-plane axis. Multiple ordinates with the same
    ``(sign_x, sign_y)`` but different ``sign_z`` share a single
    :class:`SweepDependencyGraph` instance.

    The pair :math:`(0, 0)` denotes the "pure-:math:`z`" degenerate
    set of ordinates (no streaming in the 2-D plane); the wavefront
    sweep handles those via a ``Q / Σ_t`` short-circuit and skips
    the dependency graph entirely. No graph is built for the
    :math:`(0, 0)` label.

    Frozen + slotted: hashable for use as a ``dict`` key (the SNMesh
    holds ``Mapping[OctantLabel, SweepDependencyGraph]``).
    """

    sign_x: int
    sign_y: int

    def __post_init__(self) -> None:
        if self.sign_x not in (-1, 0, +1):
            raise ValueError(f"sign_x must be in {{-1, 0, +1}}; got {self.sign_x}")
        if self.sign_y not in (-1, 0, +1):
            raise ValueError(f"sign_y must be in {{-1, 0, +1}}; got {self.sign_y}")

    @property
    def streams_in_2d(self) -> bool:
        """``False`` for the pure-z degenerate label; ``True`` otherwise."""
        return not (self.sign_x == 0 and self.sign_y == 0)


# ═══════════════════════════════════════════════════════════════════════
# SweepDependencyGraph — per-octant causal cell DAG (2-D Cartesian)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class SweepDependencyGraph:
    r"""Per-octant causal cell DAG for the 2-D Cartesian wavefront sweep.

    The DAG topology is a regular pattern: cells :math:`(i, j)` and
    :math:`(i', j')` are mutually independent at the same topological
    level iff :math:`i + j = i' + j'` (under the canonical
    :math:`(+1, +1)` orientation; reversed indices for negative-sign
    octants). This collapses to ``nx + ny - 1`` topological levels,
    each carrying an anti-diagonal of cells whose updates are
    mutually independent.

    The closed-form precompute on a Cartesian grid lives entirely
    inside :meth:`from_cartesian_2d` and never appears in the sweep
    loop. This is structural, not hand-rolled — the "library
    version" (a generic topological-sort over an explicit DAG) would
    be over-engineering for a regular pattern that collapses to ~5
    lines of arange + diagonal extraction.

    The :meth:`apply` method does NOT inline the per-cell update
    math — it dispatches to ``cell_update.update_batch(slice_args)``,
    so the WDD / LD / EC / Step closure is pluggable. This aligns
    with the strategy contract in
    :mod:`orpheus.sn.spatial.cell_update`.

    Attributes
    ----------
    label :
        :class:`OctantLabel` carried for diagnostics + key matching.
    levels :
        Tuple of ``(ii, jj)`` pairs, one per topological level. Each
        is a pair of ``int`` arrays of equal length ``n_diag``
        carrying the cell-i / cell-j indices on that level. Walked
        in topological order.
    face_in_x, face_out_x, face_in_y, face_out_y :
        Per-axis face-index offsets. ``face_in`` is the offset added
        to ``ii`` (or ``jj``) to get the incoming-face index;
        ``face_out`` for the outgoing face. For ``sx ≥ 0``:
        ``face_in_x = 0`` (the "lower" x-face), ``face_out_x = 1``
        (the "upper"). For ``sx < 0``: reversed.

    Notes
    -----
    Frozen + (no slots — frozen suffices for hashability and
    immutability; ``slots=True`` would conflict with ``levels`` being
    a tuple-of-array, which numpy can store in slot but is awkward).

    Iteration over ``levels`` is the topological order of the upwind
    DAG — structural, not a smoking gun. Within each level the cell
    update is delegated to :meth:`CellUpdateBase.update_batch`, so
    the ordinate axis (``N_oct``) is internal to every ``apply``
    call: there is no ``for n in range(N_oct)`` anywhere in this
    module.
    """

    label: OctantLabel
    levels: tuple[tuple[np.ndarray, np.ndarray], ...]
    face_in_x: int
    face_out_x: int
    face_in_y: int
    face_out_y: int

    # ── Construction ────────────────────────────────────────────────

    @classmethod
    def from_cartesian_2d(
        cls,
        *,
        nx: int,
        ny: int,
        label: OctantLabel,
    ) -> "SweepDependencyGraph":
        r"""Build the per-octant anti-diagonal schedule for an ``nx × ny`` grid.

        The cells on topological level ``k`` are those satisfying
        ``local_i + local_j == k`` under the per-octant orientation.
        For ``sign_x = +1``: ``local_i`` walks 0 → nx-1 (forward).
        For ``sign_x = -1``: ``local_i`` walks nx-1 → 0 (reversed).
        Same for ``sign_y``.

        Parameters
        ----------
        nx, ny :
            Cell counts on the Cartesian grid.
        label :
            Octant signature with ``sign_x ≠ 0 OR sign_y ≠ 0`` (the
            DAG is undefined for the pure-z degenerate label;
            building one raises ``ValueError``).

        Raises
        ------
        ValueError
            If ``label.streams_in_2d`` is ``False`` (i.e.
            ``sign_x == 0 == sign_y``).
        """
        if not label.streams_in_2d:
            raise ValueError(
                f"Cannot build a SweepDependencyGraph for the pure-z "
                f"degenerate octant {label!r}; this octant has no "
                "in-plane streaming and the wavefront sweep handles "
                "it via a Q/Σ_t short-circuit."
            )
        sx, sy = label.sign_x, label.sign_y
        # Face index offsets — match the legacy ``_diag_cache`` build at
        # ``orpheus.sn.sweep._sweep_2d_wavefront`` (lines 766-785).
        face_in_x = 0 if sx >= 0 else 1
        face_out_x = 1 if sx >= 0 else 0
        face_in_y = 0 if sy >= 0 else 1
        face_out_y = 1 if sy >= 0 else 0

        # Per-octant cell-index arrays. Forward arange for sx >= 0,
        # reversed for sx < 0 — same pattern as legacy.
        ix_arr = np.arange(nx) if sx >= 0 else np.arange(nx)[::-1]
        iy_arr = np.arange(ny) if sy >= 0 else np.arange(ny)[::-1]

        # Anti-diagonal extraction: cells where ``local_i + local_j ==
        # k`` form level k. ``ix_arr[local_i]`` and ``iy_arr[local_j]``
        # then translate to the global cell indices in the per-octant
        # traversal order.
        levels: list[tuple[np.ndarray, np.ndarray]] = []
        for k in range(nx + ny - 1):
            i_start = max(0, k - ny + 1)
            i_end = min(nx - 1, k)
            local_i = np.arange(i_start, i_end + 1)
            local_j = k - local_i
            levels.append((ix_arr[local_i], iy_arr[local_j]))

        return cls(
            label=label,
            levels=tuple(levels),
            face_in_x=face_in_x,
            face_out_x=face_out_x,
            face_in_y=face_in_y,
            face_out_y=face_out_y,
        )

    # ── Public API ─────────────────────────────────────────────────

    def apply(
        self,
        *,
        cell_update: CellUpdateBase,
        psi_x_octant: np.ndarray,        # (N_oct, ng, nx+1, ny) — mutated in place
        psi_y_octant: np.ndarray,        # (N_oct, ng, nx, ny+1) — mutated in place
        Q_octant: np.ndarray,            # (N_oct or 1, ng, nx, ny)
        sig_t: np.ndarray,               # (ng, nx, ny)
        str_x_octant: np.ndarray,        # (N_oct, nx)
        str_y_octant: np.ndarray,        # (N_oct, ny)
        weights_octant: np.ndarray,      # (N_oct,)
        angular_flux_octant: np.ndarray, # (N_oct, ng, nx, ny) — written
        scalar_flux_buf: np.ndarray,     # (ng, nx, ny) — accumulated into
    ) -> None:
        r"""Walk the topological levels and accumulate scalar + angular flux.

        Forward-substitute along the DAG, vectorised over
        ``(N_oct, n_diag, ng)`` per level. The ordinate axis
        (``N_oct``) is INTERNAL to every ``apply`` call: there is no
        ``for n in range(N_oct)`` anywhere in this method.

        Mutation contract
        -----------------

        * ``psi_x_octant`` / ``psi_y_octant`` — outgoing face fluxes
          are scattered into these buffers in place by
          :meth:`CellUpdateBase.update_batch`. Caller is responsible
          for seeding the incoming-face entries (BC apply happens
          one level up, in the wavefront sweep) and for scattering
          the post-sweep buffers back into the persistent BC state.
        * ``angular_flux_octant`` — written at every level's cells.
          Caller is responsible for scattering back into the global
          ``angular_flux`` buffer keyed by octant indices.
        * ``scalar_flux_buf`` — accumulated into via
          ``scalar_flux_buf[ii, jj, :] += weighted_sum_psi_avg``.
          Caller seeds the buffer (typically zero on the first
          octant of an outer iteration) and reads the result at the
          end of all octants.

        Parameters
        ----------
        cell_update :
            The closure strategy. Must override
            :meth:`update_batch` (DD does; Step / LD do not yet).
        psi_x_octant, psi_y_octant :
            Per-octant face buffers; mutated in place.
        Q_octant :
            Per-octant per-cell volumetric source, already
            weight-normalised. Leading axis is ``1`` for
            isotropic-only sweeps or ``N_oct`` when an aniso
            component is folded in.
        sig_t :
            Per-cell per-group total cross section.
        str_x_octant, str_y_octant :
            Per-octant streaming coefficients ``2|μ_axis|/Δaxis``.
        weights_octant :
            Per-octant ordinate weights for scalar-flux reduction.
        angular_flux_octant :
            Per-octant angular-flux output buffer.
        scalar_flux_buf :
            Global scalar-flux accumulator.
        """
        for ii, jj in self.levels:
            slice_args = SweepCellSlice(
                ii=ii, jj=jj,
                face_in_x_idx=ii + self.face_in_x,
                face_out_x_idx=ii + self.face_out_x,
                face_in_y_idx=jj + self.face_in_y,
                face_out_y_idx=jj + self.face_out_y,
                psi_x=psi_x_octant, psi_y=psi_y_octant,
                Q=Q_octant, sig_t=sig_t,
                str_x=str_x_octant, str_y=str_y_octant,
            )
            # Issue #196 PR-INDEX-5: psi_avg returned principled
            # ``(N_oct, ng, n_diag)`` (ordinate, group, anti-diagonal).
            psi_avg = cell_update.update_batch(slice_args)  # (N_oct, ng, n_diag)
            # Scalar-flux accumulation — sum over ordinates with
            # quadrature weights. ``weights_octant`` is already
            # weight-normalised by the caller (the wavefront sweep
            # divides by Σw before invoking apply), so this is a
            # plain weighted sum.
            scalar_flux_buf[:, ii, jj] += np.einsum(
                "ngd,n->gd", psi_avg, weights_octant,
            )
            # Angular-flux scatter — write all ``(N_oct, ng, n_diag)``
            # values into the per-octant principled angular-flux buffer.
            # Indices ``ii`` and ``jj`` broadcast naturally with the
            # leading ``:`` (octant) and the second ``:`` (group).
            angular_flux_octant[:, :, ii, jj] = psi_avg


__all__ = [
    "OctantLabel",
    "SweepDependencyGraph",
]

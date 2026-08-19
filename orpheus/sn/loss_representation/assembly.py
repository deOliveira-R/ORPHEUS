r"""Per-ordinate structural assembly of the Cartesian SN streaming block.

The SN arm of the assembly mode (stencil-assembly campaign 2b): emit, per
``(ordinate, group)``, the sparse matrix :math:`M` of the within-group
``L(+C)`` action over the bulk spatial DOFs — the SAME per-cell closure
coefficients the sweep (``solve``) marches and the matvec (``apply``)
evaluates, consumed a third way as ``(row, col, value)``.

What the block IS (and is not)
==============================

``M`` acts on ONE ordinate's, ONE group's bulk DOFs — the C-ravel of the
spatial axes ⊗ the closure's per-cell moment vector (DOF index
``cell·cm + moment`` with ``cm = cell_moment_count(per_axis, d)``, the
#253 single-sourced count; DD's ``cm = 1`` degenerates to the plain cell
ravel) — at **zero boundary inflow**: the Cartesian ``L + C`` is
per-ordinate per-group block-diagonal on the bulk (no angular or energy
coupling lives in it), and the trace coupling ``A_bs`` (the inflow seed)
is a separate composite block deliberately NOT emitted here — the
consumers of this block (the walk-order triangularity gates, the sweep ≡
forward-substitution discharge of #284, DSA's ``R·A·P`` moment
reduction, #2/#200) all pose it at zero inflow. Scope: **Cartesian
only** — a curvilinear ordinate couples its angular neighbor through the
Morel–Montry closure (and the #282 lagged pole seed is precisely a
walk-order back edge), so the per-ordinate block factorization does not
exist there; :meth:`SNMesh.streaming`'s own Cartesian gate enforces this
honestly.

The one-source mechanism: coefficient extraction BY KERNEL PROBING
==================================================================

The emitter carries NO stencil spelling of its own — not DD's diamond
fold, not LD's UBLD blocks. Every coefficient is extracted through the
scheme's **production matvec kernel**
(:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
— the same body ``apply`` runs, fed by the same
:meth:`SNMesh.streaming` raw per-axis data), exploiting the scheme's
declared linearity (``is_linear``): the residual and the face
reconstruction are affine in ``(ψ̄, ψ_in, Q)``,

.. math::

    r = \bar A\,\vec\psi - \sum_b C_b\,\psi^{\rm in}_b - \bar Q,
    \qquad
    \psi^{\rm out}_a = T_a\,\vec\psi + \sum_b \alpha_{ab}\,\psi^{\rm in}_b,

so unit probes per moment slot read every block exactly (unit inputs
make the kernel's products coefficient reads):

* cell probes ``ψ̄ = e_m ⊗ 1_cells`` (``cm`` of them) → ``r`` = column
  ``m`` of the per-cell diagonal block :math:`\bar A` and ``ψ_out_a`` =
  column ``m`` of the outgoing-trace map :math:`T_a`;
* face probes ``ψ_in_b = e_f ⊗ 1_cells`` (``d·fm`` of them, ``fm =
  face_moment_count``) → ``r`` = −column ``f`` of the inflow coupling
  :math:`C_b` and ``ψ_out_a`` = column ``f`` of the face-memory map
  :math:`\alpha_{ab}`.

The SAME extraction serves every batched-kernel closure: DD's scalar
faces give ``1×1`` blocks (:math:`\bar A` = the diamond denominator,
:math:`C` = the ``2g`` coupling, :math:`T = 1/w = 2`, :math:`\alpha =
-1` — the face MEMORY that makes DD's assembled block carry dense
upstream chains), LD's bilinear faces give the UBLD ``2^d × 2^d`` cell
system with the upwind trace (:math:`\alpha = 0` — no memory, block
tri-diagonal sparsity). A sign flip ANYWHERE in the shared kernel
algebra (the diamond fold, the UBLD assembly, the trace reconstruction)
therefore moves solve, apply, AND assemble together — the L16
one-source teeth bite by construction, with zero parallel spelling to
drift.

The symbolic walk
=================

With the blocks in hand, the emitter walks the octant's
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph`
in topological order carrying, per in-flight face, the face trace's ROW
BLOCK — its linear map from the bulk DOFs (dense ``(fm, n_dof)``;
boundary inflow faces are the zero block). Per cell ``c``:

.. math::

    \text{rows}(c) = \bar A_c\,E_c - \sum_b C_{b,c}\,R^{\rm in}_b,
    \qquad
    R^{\rm out}_a = T_{a,c}\,E_c + \sum_b \alpha_{ab,c}\,R^{\rm in}_b ,

where :math:`E_c` selects the cell's own DOF columns. The assembled
block is (block-)triangular in walk order — the sweep IS forward
substitution on it, which is the object-level content of #284. The
all-zero degenerate label (pure-z ordinates over a lower-D mesh) has no
graph and no face threading — its block is the cell-probe diagonal
alone (the walk's ``Q/Σ_t`` short-circuit in matrix form).

.. seealso::

   ``.claude/plans/archive/assembly_mode_crosswalk.md`` — the 2b design
   contract (DOF numbering, gate specs, the one-source table).
   ``sweep_graph.py`` "Future direction" note — this module realizes
   the per-octant assembly it anticipated (per-ordinate here; the
   octant batching is a pure vectorization the first large-scale
   consumer can add).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterator

import numpy as np
from scipy import sparse

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.moment_layout import (
    cell_moment_count,
    face_moment_count,
    face_moment_tail,
)
from orpheus.sn.loss_representation.sweep_graph import (
    OctantLabel,
    SweepDependencyGraph,
)
from orpheus.sn.loss_representation.sweep_schedule import _octant_sweep
from orpheus.transport.spatial._ubld import octant_moment_frame_signs

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh


__all__ = [
    "assemble_ordinate_blocks",
    "ordinate_walk_order",
]


def _ordinate_label(sn_mesh: "SNMesh", ordinate: int) -> OctantLabel:
    r"""The in-plane octant label of one ordinate — through the SAME
    quadrature partition + projection the sweep schedule consumes
    (:func:`~orpheus.sn.loss_representation.sweep_schedule._octant_sweep`,
    the sole in-plane projection site)."""
    for entry in sn_mesh.quad.octants:
        sweep = _octant_sweep(entry, sn_mesh.ndim)
        if ordinate in sweep.indices:
            return sweep.label
    raise IndexError(
        f"ordinate {ordinate} not found in the quadrature octant "
        f"partition (N = {sn_mesh.quad.n_ordinates})."
    )


def _iter_cells_in_walk_order(
    graph: SweepDependencyGraph,
) -> Iterator[tuple[tuple[int, ...], int]]:
    r"""Yield ``(cell multi-index, flat C-ravel index)`` in topological order.

    The SINGLE iteration spelling shared by the emitter and
    :func:`ordinate_walk_order` — the permutation the triangularity
    gates apply IS the order the emitter walked, by construction.
    """
    shape = graph.spatial_shape
    for level in graph.levels:
        for cell in zip(*(idx.tolist() for idx in level)):
            yield cell, int(np.ravel_multi_index(cell, shape))


def ordinate_walk_order(sn_mesh: "SNMesh", ordinate: int) -> np.ndarray:
    r"""Flat bulk CELL indices in the ordinate's sweep (topological) order.

    The cell-level permutation of the walk-order triangularity claim:
    for a streaming Cartesian ordinate the assembled block is
    (block-)lower-triangular under it — the object-level form of "the
    sweep is a forward substitution" (#284). For a multi-moment closure
    the DOF-level permutation is the block expansion
    ``concat(c·cm + arange(cm) for c in order)`` (the DOF layout is
    cell-major, moment-minor — the bulk C-ravel). The degenerate
    all-zero label (no streaming) returns the identity order (its block
    is cell-diagonal — trivially triangular in any order).
    """
    label = _ordinate_label(sn_mesh, ordinate)
    n_cells = int(np.prod(sn_mesh.spatial_shape))
    if not label.streams:
        return np.arange(n_cells)
    graph = SweepDependencyGraph.for_shape(sn_mesh.spatial_shape)[label]
    return np.fromiter(
        (flat for _, flat in _iter_cells_in_walk_order(graph)),
        dtype=np.intp,
        count=n_cells,
    )


def _probe_coefficient_blocks(
    scheme,
    s_axes_k: tuple[np.ndarray, ...],
    reaction_flat: np.ndarray,
    *,
    cm: int,
    fm: int,
    n_cells: int,
):
    r"""Extract the per-cell coefficient blocks by unit kernel probes.

    Returns ``(diag, inflow, trace, memory)`` with shapes

    * ``diag``   ``(ng, n_cells, cm, cm)`` — :math:`\bar A` per cell;
    * ``inflow`` ``d × (ng, n_cells, cm, fm)`` — :math:`C_b`;
    * ``trace``  ``d × (ng, n_cells, fm, cm)`` — :math:`T_a`;
    * ``memory`` ``d×d × (ng, n_cells, fm, fm)`` — :math:`\alpha_{ab}`
      (``memory[a][b]``; DD carries the diamond ``−1`` on ``a == b``,
      LD's upwind trace carries all zeros).

    Shape conventions follow the batched-kernel contract exactly: the
    cell moment axis is PRESENT iff ``cm > 1`` and the face moment tail
    iff ``fm > 1`` (:func:`face_moment_tail` — the #240 D5b cochain
    policy the walks share).
    """
    d = len(s_axes_k)
    ng = reaction_flat.shape[0]
    zeros_q = np.zeros((1, ng, n_cells))
    cell_tail = (cm,) if cm > 1 else ()
    face_tail = face_moment_tail(fm)

    def cell_field(slot: "int | None") -> np.ndarray:
        field = np.zeros((1, ng, n_cells) + cell_tail)
        if slot is not None:
            if cell_tail:
                field[..., slot] = 1.0
            else:
                field[...] = 1.0
        return field

    def face_fields(axis: "int | None", slot: int = 0):
        fields = []
        for b in range(d):
            field = np.zeros((1, ng, n_cells) + face_tail)
            if b == axis:
                if face_tail:
                    field[..., slot] = 1.0
                else:
                    field[...] = 1.0
            fields.append(field)
        return tuple(fields)

    def as_cell_columns(residual: np.ndarray) -> np.ndarray:
        # (1, ng, n_cells[, cm]) → (ng, n_cells, cm)
        r = residual[0]
        return r[..., None] if not cell_tail else r

    def as_face_columns(face_out: np.ndarray) -> np.ndarray:
        # (1, ng, n_cells[, fm]) → (ng, n_cells, fm)
        f = face_out[0]
        return f[..., None] if not face_tail else f

    diag = np.empty((ng, n_cells, cm, cm))
    trace = [np.empty((ng, n_cells, fm, cm)) for _ in range(d)]
    for m in range(cm):
        residual, faces_out = scheme.residual_kernel_batch(
            psi_bar=cell_field(m), psi_in=face_fields(None),
            s_axes=s_axes_k, reaction_xs=reaction_flat, Q_cells=zeros_q,
        )
        diag[:, :, :, m] = as_cell_columns(residual)
        for a in range(d):
            trace[a][:, :, :, m] = as_face_columns(faces_out[a])

    inflow = [np.empty((ng, n_cells, cm, fm)) for _ in range(d)]
    memory = [
        [np.empty((ng, n_cells, fm, fm)) for _ in range(d)] for _ in range(d)
    ]
    for b in range(d):
        for f in range(fm):
            residual, faces_out = scheme.residual_kernel_batch(
                psi_bar=cell_field(None), psi_in=face_fields(b, f),
                s_axes=s_axes_k, reaction_xs=reaction_flat, Q_cells=zeros_q,
            )
            inflow[b][:, :, :, f] = -as_cell_columns(residual)
            for a in range(d):
                memory[a][b][:, :, :, f] = as_face_columns(faces_out[a])
    return diag, inflow, trace, memory


def assemble_ordinate_blocks(
    sn_mesh: "SNMesh",
    ordinate: int,
    *,
    include_collision: bool = True,
) -> tuple[SparseAssembledOperator, ...]:
    r"""Emit the per-group sparse blocks of ``L(+C)`` for one ordinate.

    Parameters
    ----------
    sn_mesh :
        A CARTESIAN SN phase space (slab or 2-D; the
        :meth:`SNMesh.streaming` accessor is the Cartesian gate) with a
        LINEAR batched-kernel scheme (``is_linear`` — the extraction
        precondition; DD and LD both declare it). The scheme,
        quadrature, widths, and cross sections are all read from the
        mesh — the same sources the sweep walk consumes.
    ordinate :
        Index into the quadrature's ordinate axis.
    include_collision :
        ``True`` (default): the sweep-invertible ``L + C`` block (the
        object ``(L+C).solve`` inverts — the probes' ``reaction_xs`` is
        the mesh's :math:`\Sigma_t` field). ``False``: the pure
        streaming block (``reaction_xs = 0``; NOTE it is singular for a
        non-streaming degenerate ordinate, honestly: an all-zero
        matrix).

    Returns
    -------
    tuple[SparseAssembledOperator, ...]
        One ``(n_cells·cm × n_cells·cm)`` block per energy group, over
        the bulk DOF layout ``cell·cm + moment`` (the C-ravel of
        ``(*spatial, cm)``), at zero boundary inflow (module
        docstring). ``block.apply(x)`` reproduces the production
        matvec's ``(ordinate, group)`` bulk row on a bulk field that is
        zero outside that row.
    """
    spatial_shape = tuple(int(n) for n in sn_mesh.spatial_shape)
    n_cells = int(np.prod(spatial_shape))
    d = len(spatial_shape)
    scheme = sn_mesh.scheme
    if not type(scheme).is_linear:
        from orpheus.numerics.operator import MissingAssembly

        raise MissingAssembly(
            f"assemble_ordinate_blocks: coefficient extraction by unit "
            f"probes requires a LINEAR closure (is_linear); "
            f"{type(scheme).__name__} declares otherwise."
        )
    per_axis = int(scheme.spatial_basis_per_axis)
    cm = cell_moment_count(per_axis, d)
    fm = face_moment_count(per_axis, d)
    n_dof = n_cells * cm

    # ── The raw data every mode shares ────────────────────────────────
    sigma_t = np.asarray(
        sn_mesh.material_xs_field().total_cross_section_field.values, float,
    )                                                   # (ng, *spatial)
    ng = sigma_t.shape[0]
    reaction_flat = (
        sigma_t.reshape(ng, n_cells)
        if include_collision
        else np.zeros((ng, n_cells))
    )
    # Per-axis raw down-face streaming for THIS ordinate, gathered per
    # flat cell: streaming(a) is the production (N, n_a) accessor the
    # DAG walk reads as str_axes[a] (#240: the scheme applies its own
    # closure factor, never this data).
    cell_index = np.indices(spatial_shape).reshape(d, n_cells)
    s_axes_k = tuple(
        np.asarray(sn_mesh.streaming(a), float)[ordinate][cell_index[a]][
            None, None, :
        ]
        for a in range(d)
    )

    diag, inflow, trace, memory = _probe_coefficient_blocks(
        scheme, s_axes_k, reaction_flat, cm=cm, fm=fm, n_cells=n_cells,
    )

    # ── The symbolic walk (or the degenerate cell-diagonal) ───────────
    label = _ordinate_label(sn_mesh, ordinate)

    # The kernel produces/consumes moments in the per-ordinate SWEEP
    # frame; the bulk field lives in the GLOBAL frame, and the walks
    # convert with the ``_reframe`` involution. The assembled block must
    # be the GLOBAL-frame operator, so conjugate by the SAME sign vector
    # the walk uses (one source: ``octant_moment_frame_signs`` — ``None``
    # for a single-moment closure, where the frames coincide):
    # ``M_global = Φ · M_sweep · Φ`` with ``Φ = diag(signs per cell)``.
    frame_signs = octant_moment_frame_signs(label.signs, per_axis)
    dof_signs = (
        None if frame_signs is None else np.tile(frame_signs, n_cells)
    )

    def _to_global_frame(
        vals: np.ndarray, rows: np.ndarray, cols: np.ndarray,
    ) -> np.ndarray:
        if dof_signs is None:
            return vals
        return vals * dof_signs[rows] * dof_signs[cols]
    if not label.streams:
        # The Q/Σ_t short-circuit's operator: the cell-diagonal blocks
        # alone (exactly the probe diagonal at zero couplings).
        blocks = []
        base = np.arange(n_cells)[:, None, None] * cm
        rows = (base + np.arange(cm)[:, None]).ravel()
        cols = (base + np.arange(cm)[None, :]).ravel()
        for g in range(ng):
            vals = _to_global_frame(diag[g].ravel(), rows, cols)
            blocks.append(
                SparseAssembledOperator(
                    sparse.coo_array(
                        (vals, (rows.copy(), cols.copy())),
                        shape=(n_dof, n_dof),
                    )
                )
            )
        return tuple(blocks)

    graph = SweepDependencyGraph.for_shape(spatial_shape)[label]
    blocks = []
    for g in range(ng):
        # In-flight face row blocks, keyed by (axis, face multi-index):
        # each is the face trace's dense (fm, n_dof) map. Boundary
        # inflow faces are absent ⟹ the zero block (zero-inflow posing).
        face_blocks: dict[tuple[int, tuple[int, ...]], np.ndarray] = {}
        all_rows: list[np.ndarray] = []
        all_cols: list[np.ndarray] = []
        all_vals: list[np.ndarray] = []
        for cell, flat in _iter_cells_in_walk_order(graph):
            base = flat * cm
            incoming: list["np.ndarray | None"] = []
            cell_rows = np.zeros((cm, n_dof))
            cell_rows[:, base:base + cm] = diag[g, flat]
            for b in range(d):
                key_in = (
                    b,
                    cell[:b] + (cell[b] + graph.face_in[b],) + cell[b + 1:],
                )
                row_in = face_blocks.pop(key_in, None)
                incoming.append(row_in)
                if row_in is not None:
                    cell_rows -= inflow[b][g, flat] @ row_in
            for a in range(d):
                key_out = (
                    a,
                    cell[:a] + (cell[a] + graph.face_out[a],) + cell[a + 1:],
                )
                row_out = np.zeros((fm, n_dof))
                row_out[:, base:base + cm] = trace[a][g, flat]
                for b, row_in in enumerate(incoming):
                    if row_in is not None:
                        row_out += memory[a][b][g, flat] @ row_in
                face_blocks[key_out] = row_out
            local_rows, nonzero_cols = np.nonzero(cell_rows)
            all_rows.append(base + local_rows)
            all_cols.append(nonzero_cols)
            all_vals.append(cell_rows[local_rows, nonzero_cols])
        rows = np.concatenate(all_rows)
        cols = np.concatenate(all_cols)
        vals = _to_global_frame(np.concatenate(all_vals), rows, cols)
        blocks.append(
            SparseAssembledOperator(
                sparse.coo_array(
                    (vals, (rows, cols)), shape=(n_dof, n_dof),
                )
            )
        )
    return tuple(blocks)

r"""#310 C3 — the multi-D reverse ``walk_full`` gates (spec §5, R2a slice).

The multi-D Cartesian adjoint matvec ``(L+C)ᵀφ`` exists as the ORACLE arm:
:meth:`FullFieldWavefront.loss_action_transpose` routes through the shared
:class:`_OctantWalk` apply-transpose frame — the UNCHANGED
:meth:`SweepDependencyGraph.walk_full` over each octant's MIRROR graph
(``−signs_eff``: reversed levels + swapped face roles, the discrete face of
"the adjoint streams along −Ω") × the :class:`_CellResidualTranspose` level
operation, bottoming in the SAME scheme kernel VJP as the 1-D reverse arms
(``residual_kernel_batch_transpose``, #310 C2).

This file is the multi-D sibling of ``test_one_dim_loop_walk.py`` (spec
§5.2) plus the §5.3 structurally-independent object oracles:

* **runtime spy** — both matvec orientations exercise the ONE
  ``walk_full`` frame; direction is the level-op OBJECT + mirror-label
  DATA, observed at runtime;
* **AST tripwire** — no boolean orientation flags in the shared frames;
* **dense-``Mᵀ`` column probe** (the NEW 2-D oracle artifact) — ``M``
  probed off the FORWARD apply over the full composite basis, the reverse
  probed the same way, ``M_rev == M_fwdᵀ`` pinned as an OBJECT (the Mode-12
  stabiliser escape: a matrix equality, outside every spectral invariance
  group);
* **Euclidean pairing identity** — ``⟨Fx, w⟩ = ⟨x, Fᵀw⟩`` over the full
  composite (bulk ⊕ trace; covers the boundary defect/identity algebra);
* **d=1 cross-realization** — the SAME slab operator through TWO
  independent reverse realizations (the C2 reverse legs vs the C3
  mirror-DAG walk) is BIT-identical (same kernel, same elementwise ops,
  batching-order-free);
* **assembled-``Mᵀ``** — the CSR ``M.T @ x`` of the forward-probed
  per-ordinate bulk blocks (LAPACK-side artifact, structurally independent
  of the walk);
* **axis equivariance + mutation teeth** (M-R2-ADDRESSING + the
  M-R2-AXISSWAP partial-swap tooth) — the committed value teeth of spec
  §11, plus a MEASURED design finding each way: (a) M-R2-LEVELORDER is
  *unrepresentable* in the mirror-graph realization — the level order and
  the face roles are ONE graph object, so "reverse the addressing but not
  the levels" cannot be spelled; (b) the TOTAL axis conjugation is an
  exact no-op on het σ (2e-16) — the reverse interior is genuinely
  d-generic (no hidden per-axis code), so the representable axis-swap bug
  class is the PARTIAL swap (one tuple crossed, the Mode-2 variable
  swap), which reds the pairing O(1) where it types (square) and cannot
  even shape-check on the rectangular primary configs (L16 — the
  shape-guard is the reason nx≠ny is mandated);
* **loud deferrals** — LD-2D reverse (typed → #310 C5), the wavefront
  ``sweep_transpose`` (G-S reverse-solve, out of scope R7), and the
  Pattern-4 moment-tail backstop.

``-O``-safe (vv Mode 8): ``pytest.fail`` / ``np.testing`` only.
``foundation`` — software/algebra invariants (no theory ``:label:``).
"""
from __future__ import annotations

import ast
import inspect
import textwrap

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation import (
    CumprodScan,
    FullFieldWavefront,
    _OctantWalk,
    _reverse_octant_traversal,
)
from orpheus.sn.loss_representation.assembly import assemble_ordinate_blocks
from orpheus.sn.loss_representation.sweep_graph import (
    SweepDependencyGraph,
    _CellResidualTranspose,
)
from orpheus.sn.loss_representation import _ApplyOperands
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from tests.sn.operators.test_g_adjoint_reciprocity import (
    _make_ld_slab,
    _make_slab,
    _random_composite,
)

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Builders — het σ + non-uniform h + rectangular nx≠ny (spec §5.5 / §9)
# ═══════════════════════════════════════════════════════════════════════


def _cart2d_probe_mesh() -> SNMesh:
    """Small rectangular (nx=3 ≠ ny=2) NON-UNIFORM 2-material cart2d mesh —
    the dense-probe config.  Rectangular + non-uniform h + het σ makes
    ``A ≠ Aᵀ`` observable (L16: a square-uniform-symmetric config is
    transpose-BLIND to axis-swap/DOF-transposition classes)."""
    geom = Mesh2D(
        edges_x=np.array([0.0, 0.4, 1.1, 2.0]),
        edges_y=np.array([0.0, 0.7, 1.5]),
        mat_map=np.array([[0, 1], [1, 0], [0, 0]]),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    return SNMesh(
        geom, Quadrature.level_symmetric(2),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )


def _cart2d_square_uniform_mesh() -> SNMesh:
    """Square (3×3) UNIFORM mesh for the axis-conjugation tooth — the config
    on which the x↔y swap is a symmetry of the GEOMETRY (the quadrature's S2
    ordinates have ``μ_x = μ_y`` exactly), isolating what the mutation
    breaks: the octant-identity coupling, not the mesh."""
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, 4),
        edges_y=np.linspace(0.0, 2.0, 4),
        mat_map=np.zeros((3, 3), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    return SNMesh(
        geom, Quadrature.level_symmetric(2), {0: get_mixture("A", "2g")},
    )


def _het_sigma(sn: SNMesh, rng: np.random.Generator) -> np.ndarray:
    """Heterogeneous (space × group) positive σ_t for the rep-level calls."""
    return 0.4 + rng.random((2, *sn.spatial_shape))


# ═══════════════════════════════════════════════════════════════════════
# Composite ⇄ flat codec (local, used consistently on BOTH directions —
# the transpose claim is codec-invariant under any fixed DOF bijection)
# ═══════════════════════════════════════════════════════════════════════


def _zero_composite(sn: SNMesh) -> FullField:
    return FullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
    )


def _flatten(field, faces: tuple[str, ...]) -> np.ndarray:
    parts = [np.asarray(field.interior.values).ravel()]
    parts += [np.asarray(field.boundary.face_view(f)).ravel() for f in faces]
    return np.concatenate(parts)


def _basis_size(sn: SNMesh, faces: tuple[str, ...]) -> int:
    z = _zero_composite(sn)
    n = int(np.asarray(z.interior.values).size)
    for f in faces:
        n += int(np.asarray(z.boundary.face_view(f)).size)
    return n


def _basis_composites(sn: SNMesh, faces: tuple[str, ...]):
    """Unit composites in the SAME DOF order ``_flatten`` reads."""
    z = _zero_composite(sn)
    bulk_shape = np.asarray(z.interior.values).shape
    for idx in np.ndindex(*bulk_shape):
        e = _zero_composite(sn)
        e.interior.values[idx] = 1.0
        yield e
    for f in faces:
        view_shape = np.asarray(z.boundary.face_view(f)).shape
        for idx in np.ndindex(*view_shape):
            e = _zero_composite(sn)
            e.boundary.face_view(f)[idx] = 1.0
            yield e


def _pairing(a, b, faces: tuple[str, ...]) -> float:
    tot = float(np.sum(np.asarray(a.interior.values) * np.asarray(b.interior.values)))
    for f in faces:
        tot += float(np.sum(
            np.asarray(a.boundary.face_view(f)) * np.asarray(b.boundary.face_view(f))
        ))
    return tot


def _pairing_defect(sn: SNMesh, rep, sig: np.ndarray, rng) -> float:
    """Relative defect of ``⟨Fx, w⟩ − ⟨x, Fᵀw⟩`` on random full composites."""
    faces = tuple(sn.angular_trace.face_names)
    x = _random_composite(sn, rng)
    w = _random_composite(sn, rng)
    lhs = _pairing(rep.loss_action(sig, x), w, faces)
    rhs = _pairing(x, rep.loss_action_transpose(sig, w), faces)
    return abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)


# ═══════════════════════════════════════════════════════════════════════
# §5.2 — the runtime spy + the AST tripwire (the one-walk claim)
# ═══════════════════════════════════════════════════════════════════════


def test_reverse_matvec_routes_through_the_shared_frame(monkeypatch):
    """[L0 structural] the multi-D reverse routes through the SAME
    ``walk_full`` + ``_interior_walk`` frames the forward uses, with the
    ``_CellResidualTranspose`` level op — and the FORWARD control leg never
    touches the transpose op (the two-direction one-walk matrix, Mode-11
    wrap sentinel: an in-process counter, not a green-only claim)."""
    walk_hits: list[str] = []
    op_hits: list[str] = []
    frame_hits: list[str] = []
    real_walk = SweepDependencyGraph.walk_full
    real_cell = _CellResidualTranspose.cell
    real_frame = _OctantWalk._interior_walk

    def walk_spy(self, **kwargs):
        walk_hits.append("walked")
        return real_walk(self, **kwargs)

    def cell_spy(self, cell_idx, **kwargs):
        op_hits.append("cell")
        return real_cell(self, cell_idx, **kwargs)

    def frame_spy(self, sweeps, **kwargs):
        frame_hits.append("frame")
        return real_frame(self, sweeps, **kwargs)

    monkeypatch.setattr(SweepDependencyGraph, "walk_full", walk_spy)
    monkeypatch.setattr(_CellResidualTranspose, "cell", cell_spy)
    monkeypatch.setattr(_OctantWalk, "_interior_walk", frame_spy)

    rng = np.random.default_rng(20260724)
    sn = _cart2d_probe_mesh()
    sig = _het_sigma(sn, rng)
    rep = FullFieldWavefront(sn)
    phi = _random_composite(sn, rng)

    # (1) the REVERSE matvec: all three sentinels fire.
    _ = rep.loss_action_transpose(sig, phi)
    if not (walk_hits and op_hits and frame_hits):
        pytest.fail(
            "the reverse matvec did NOT route through the shared frames "
            f"(walk_full hits={len(walk_hits)}, "
            f"_CellResidualTranspose.cell hits={len(op_hits)}, "
            f"_interior_walk hits={len(frame_hits)}) — the C3 carve "
            "regressed into a private reverse walk."
        )

    # (2) the FORWARD control: the same walk fires, the transpose op NEVER.
    walk_hits.clear(); op_hits.clear(); frame_hits.clear()
    _ = rep.loss_action(sig, phi)
    if not (walk_hits and frame_hits):
        pytest.fail(
            "the forward matvec left the shared frames — the one-walk "
            "claim broke on the forward side."
        )
    if op_hits:
        pytest.fail(
            f"_CellResidualTranspose.cell fired {len(op_hits)}× during the "
            "FORWARD matvec — the direction objects leaked across "
            "orientations."
        )


def test_reverse_walk_is_orientation_object_not_boolean():
    """[L0 structural] the shared multi-D frames fork on OBJECTS (mirror
    labels + level-op) — never a boolean orientation flag (the
    ``test_one_dim_loop_walk`` rule's multi-D sibling; spec §5.2 M-R2-SPY).
    Source-inspection via AST identifiers, so docstrings NAMING the
    anti-pattern don't trip it.  ``-O``-safe."""
    smells = {
        "is_solve", "is_apply", "is_matvec",
        "is_adjoint", "is_forward", "is_transpose", "is_reverse",
    }
    offenders: set[str] = set()
    for obj in (
        SweepDependencyGraph,
        _CellResidualTranspose,
        _OctantWalk,
        _reverse_octant_traversal,
    ):
        tree = ast.parse(textwrap.dedent(inspect.getsource(obj)))
        identifiers = {
            node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
        } | {
            node.attr for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)
        } | {
            node.arg for node in ast.walk(tree)
            if isinstance(node, (ast.arg, ast.keyword))
            and node.arg is not None
        }
        offenders |= identifiers & smells
    if offenders:
        pytest.fail(
            f"the multi-D walk frames carry boolean orientation flag(s) "
            f"{sorted(offenders)} — orientation MUST be carried by the "
            "mirror-label data + the level-op objects (coding-elegance "
            "Smell #3)."
        )


# ═══════════════════════════════════════════════════════════════════════
# §5.3 — the structurally-independent object oracles
# ═══════════════════════════════════════════════════════════════════════


def test_dense_mt_2d_column_probe_pins_the_object():
    """[L0 object] the NEW 2-D dense-``Mᵀ`` oracle: ``M`` column-probed off
    the FORWARD apply over the FULL composite basis (bulk ⊕ trace), the
    reverse probed identically — ``M_rev == M_fwdᵀ`` as a MATRIX equality.

    Mode-12: this pins the OBJECT, outside every spectral invariance group
    (``eig(Mᵀ) = eig(M)`` makes any spectral functional designed-green on
    the whole transpose mutation class).  Config per spec §5.5: rectangular
    nx≠ny + non-uniform h + het σ, so ``M ≠ Mᵀ`` and addressing bugs are
    observable.  The full composite is mandatory — the domain-boundary
    in↔out swap lives in the trace rows (a bulk-only probe is blind to it).
    """
    rng = np.random.default_rng(20260726)
    sn = _cart2d_probe_mesh()
    sig = _het_sigma(sn, rng)
    rep = FullFieldWavefront(sn)
    faces = tuple(sn.angular_trace.face_names)
    n = _basis_size(sn, faces)

    M_fwd = np.empty((n, n))
    M_rev = np.empty((n, n))
    for k, e in enumerate(_basis_composites(sn, faces)):
        M_fwd[:, k] = _flatten(rep.loss_action(sig, e), faces)
    for k, e in enumerate(_basis_composites(sn, faces)):
        M_rev[:, k] = _flatten(rep.loss_action_transpose(sig, e), faces)

    scale = float(np.max(np.abs(M_fwd)))
    np.testing.assert_allclose(
        M_rev, M_fwd.T, rtol=1e-12, atol=1e-13 * scale,
        err_msg=(
            "the reverse walk is NOT the transpose of the forward walk "
            "(dense column-probe object mismatch)"
        ),
    )
    # The config genuinely discriminates: the operator must not be
    # accidentally symmetric (a symmetric M would null the whole gate).
    if np.allclose(M_fwd, M_fwd.T, rtol=1e-6, atol=1e-9 * scale):
        pytest.fail(
            "probe config produced a symmetric M — the dense-Mᵀ gate is "
            "vacuous on this config (L16: pick het/non-uniform/rectangular)"
        )


def test_pairing_identity_full_composite():
    """[L0 object] ``⟨Fx, w⟩ = ⟨x, Fᵀw⟩`` at machine precision on random
    full composites — the whole-surface Euclidean pairing (bulk residual +
    chain + boundary defect/identity algebra in one identity), on both the
    vacuum probe mesh and the reflective nonsquare helper mesh."""
    from tests.sn._test_helpers import cart2d_2g_nonsquare

    rng = np.random.default_rng(20260727)
    for name, sn in (
        ("vacuum 3x2 het", _cart2d_probe_mesh()),
        ("reflective 5x7", cart2d_2g_nonsquare()),
    ):
        sig = _het_sigma(sn, rng)
        rel = _pairing_defect(sn, FullFieldWavefront(sn), sig, rng)
        if rel > 1e-12:
            pytest.fail(
                f"[{name}] Euclidean pairing identity broke: rel defect "
                f"{rel:.3e} (exact-transpose claim is machine-precision)"
            )


def test_d1_cross_realization_bit_identical():
    """[L0 cross] the SAME slab operator through TWO independent reverse
    realizations — the C2 reverse legs (``CumprodScan``) vs the C3
    mirror-DAG ``walk_full`` (``FullFieldWavefront``) — is BIT-identical,
    DD and LD.

    Byte-equality is structural, not luck: both bottom in the same
    ``residual_kernel_batch_transpose`` with the same per-cell operand
    values, and every op is elementwise over the ordinate axis (the leg
    batching vs octant batching difference cannot move a bit).  If this
    ever drifts to ULP, a kernel gained a cross-ordinate reduction —
    re-derive before relaxing."""
    rng = np.random.default_rng(20260728)
    for name, (sn, sig) in (
        ("DD slab", _make_slab(ng=2)),
        ("LD slab het non-uniform", _make_ld_slab(ng=2)),
    ):
        x = _random_composite(sn, rng)
        z_scan = CumprodScan(sn).loss_action_transpose(sig, x)
        z_dag = FullFieldWavefront(sn).loss_action_transpose(sig, x)
        np.testing.assert_array_equal(
            np.asarray(z_dag.interior.values),
            np.asarray(z_scan.interior.values),
            err_msg=f"[{name}] bulk cotangents differ between realizations",
        )
        for f in sn.angular_trace.face_names:
            np.testing.assert_array_equal(
                np.asarray(z_dag.boundary.face_view(f)),
                np.asarray(z_scan.boundary.face_view(f)),
                err_msg=f"[{name}] trace cotangent differs on face {f}",
            )


def test_assembled_mt_2d_per_ordinate_block():
    """[L0 object] assembled-``Mᵀ`` (2-D): the CSR ``M.T @ x`` of each
    forward-probed per-ordinate bulk block equals the reverse walk's bulk
    output on a bulk-impulse cotangent — and the transpose stays exactly
    per-(ordinate, group) block-diagonal (no off-block leak).

    Structurally independent of the walk under test: the blocks are
    assembled by forward-kernel unit probes and transposed by scipy CSR.
    σ is the mesh's own material field — the SAME source the assembly
    reads."""
    sn = SNMesh(
        Mesh2D(
            edges_x=np.array([0.0, 0.4, 1.1, 2.1, 3.0]),
            edges_y=np.array([0.0, 0.7, 1.5, 2.0]),
            mat_map=np.array([[0, 1, 1], [1, 0, 0], [0, 0, 1], [1, 1, 0]]),
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
        ),
        Quadrature.level_symmetric(sn_order=4),
        {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
    )
    sigma = np.asarray(
        sn.material_xs_field().total_cross_section_field.values, float,
    )
    rep = FullFieldWavefront(sn)
    rng = np.random.default_rng(20260729)
    n_cells = int(np.prod(sn.spatial_shape))
    for n in range(sn.quad.n_ordinates):
        blocks = assemble_ordinate_blocks(sn, n)
        for g in range(2):
            r = rng.standard_normal(n_cells)
            w = _zero_composite(sn)
            w.interior.values[n, g] = r.reshape(sn.spatial_shape)
            z = rep.loss_action_transpose(sigma, w)
            bulk = np.asarray(z.interior.values)
            np.testing.assert_allclose(
                bulk[n, g].ravel(), blocks[g].apply_transpose(r),
                rtol=1e-12, atol=1e-13,
                err_msg=f"assembled-Mᵀ broke at ordinate {n}, group {g}",
            )
            rest = bulk.copy()
            rest[n, g] = 0.0
            np.testing.assert_array_equal(
                rest, 0.0,
                err_msg=(
                    f"reverse walk leaked off-block at (n={n}, g={g}) — "
                    "the transpose must stay per-ordinate block-diagonal"
                ),
            )


# ═══════════════════════════════════════════════════════════════════════
# §11 — the committed mutation teeth (M-R2-ADDRESSING + the axis swap)
# ═══════════════════════════════════════════════════════════════════════


def test_addressing_mutation_reds(monkeypatch):
    """[Mode-10 tooth] M-R2-ADDRESSING: dropping the mirror (forward labels
    in the reverse traversal — gather at ``face_in``, seed at the physical
    in-edge, i.e. the whole transposed addressing gone) reds the pairing
    identity O(1).  Fires under ``-O``; reverted by monkeypatch."""
    import orpheus.sn.loss_representation as lr

    rng = np.random.default_rng(20260730)
    sn = _cart2d_probe_mesh()
    sig = _het_sigma(sn, rng)

    monkeypatch.setattr(lr, "_reverse_octant_traversal", lambda sweeps: sweeps)
    rel = _pairing_defect(sn, FullFieldWavefront(sn), sig, rng)
    if rel < 1e-3:
        pytest.fail(
            f"un-mirrored reverse traversal moved the pairing defect only "
            f"{rel:.3e} — the addressing tooth has no bite (the gate would "
            "miss a forward-addressed reverse walk)"
        )


def test_reverse_interior_is_axis_equivariant(monkeypatch):
    """[L0 invariance] the reverse interior is AXIS-EQUIVARIANT: the total
    x↔y conjugation (crossed labels, face tuples, streaming axes, spatial
    axes — inputs AND outputs) is an exact no-op, on HET σ (the
    discriminating config: any hidden per-axis special case would break
    equivariance there).

    This is the symmetry-group pin of d-genericity — the reverse path
    contains no hand-rolled x/y code, so conjugating its whole world by
    the axis permutation reproduces it bit-for-machine-precision.  It is
    NOT a mutation tooth (a total conjugation is a strict symmetry of a
    d-generic FUNCTION, any config); the representable bug class is the
    PARTIAL swap — the companion tooth below."""
    rng = np.random.default_rng(20260731)
    sn = _cart2d_square_uniform_mesh()

    orig = FullFieldWavefront._loss_action_transpose_interior

    def _swap_spatial(a: np.ndarray) -> np.ndarray:
        return np.swapaxes(a, -2, -1)

    def conjugated_interior(self, operands, oct_idx, signs_addr, out_bars):
        ops2 = _ApplyOperands(
            probe=_swap_spatial(operands.probe),
            sig_t=_swap_spatial(operands.sig_t),
            str_axes=operands.str_axes[::-1],
            Q_zero=_swap_spatial(operands.Q_zero),
        )
        psi_cot, capture = orig(
            self, ops2, oct_idx, signs_addr[::-1], out_bars[::-1],
        )
        return _swap_spatial(psi_cot), capture[::-1]

    monkeypatch.setattr(
        FullFieldWavefront, "_loss_action_transpose_interior",
        conjugated_interior,
    )
    sig_het = 0.4 + rng.random((2, 3, 3))       # transpose-ASYMMETRIC
    rel = _pairing_defect(sn, FullFieldWavefront(sn), sig_het, rng)
    if rel > 1e-12:
        pytest.fail(
            f"total axis conjugation moved the reverse by {rel:.3e} on het "
            "σ — the reverse interior gained axis-specific (non-d-generic) "
            "code"
        )


def test_axis_swap_partial_mutation_reds(monkeypatch):
    """[Mode-10 tooth] M-R2-AXISSWAP as the representable bug: the PARTIAL
    swap — crossing ONE per-axis tuple (the out-face cotangents) against
    uncrossed addressing, the Mode-2 variable swap — reds the pairing
    O(1).  Demonstrated on the SQUARE mesh because that is where it TYPES:
    on the rectangular primary configs (dense-Mᵀ, assembled-Mᵀ) the same
    mutation cannot even shape-check — the L16 reason nx≠ny is mandated
    (the whole silent-value-bug class becomes a loud shape error)."""
    rng = np.random.default_rng(20260732)
    sn = _cart2d_square_uniform_mesh()
    sig = 0.4 + rng.random((2, 3, 3))

    orig = FullFieldWavefront._loss_action_transpose_interior

    def crossed_faces_interior(self, operands, oct_idx, signs_addr, out_bars):
        return orig(self, operands, oct_idx, signs_addr, out_bars[::-1])

    monkeypatch.setattr(
        FullFieldWavefront, "_loss_action_transpose_interior",
        crossed_faces_interior,
    )
    rel = _pairing_defect(sn, FullFieldWavefront(sn), sig, rng)
    if rel < 1e-3:
        pytest.fail(
            f"crossed out-face cotangent tuple moved the pairing defect "
            f"only {rel:.3e} — the partial axis-swap tooth has no bite"
        )


# ═══════════════════════════════════════════════════════════════════════
# Loud deferrals (spec §12.2 — out-of-scope stays typed and RED-loud)
# ═══════════════════════════════════════════════════════════════════════


def test_ld_2d_reverse_is_a_typed_deferral():
    """[deferral pin] the LD-2D reverse would RUN (the kernel VJP is
    d-generic and LD registers it) but is UNGATED until #310 C5 — the frame
    refuses loudly rather than emitting an unverified answer."""
    from orpheus.transport.spatial.linear_discontinuous import (
        LinearDiscontinuous,
    )

    geom = Mesh2D(
        edges_x=np.array([0.0, 0.5, 1.3, 2.0]),
        edges_y=np.array([0.0, 0.9, 2.0]),
        mat_map=np.zeros((3, 2), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    sn = SNMesh(
        geom, Quadrature.level_symmetric(2),
        {0: get_mixture("A", "2g")},
        scheme=LinearDiscontinuous(),
    )
    rng = np.random.default_rng(20260801)
    phi = _random_composite(sn, rng)
    with pytest.raises(NotImplementedError, match="#310 C5"):
        FullFieldWavefront(sn).loss_action_transpose(_het_sigma(sn, rng), phi)


def test_wavefront_solve_transpose_still_raises():
    """[deferral pin] the multi-D ``sweep_transpose`` (the G-S reverse-SOLVE,
    out-of-scope R7) stays a typed raise — C3 lands the matvec transpose
    ONLY, and must not silently un-defer the solve arm."""
    sn = _cart2d_probe_mesh()
    with pytest.raises(NotImplementedError, match="reverse-scan"):
        FullFieldWavefront(sn).sweep_transpose(
            np.zeros((sn.quad.N, 2, *sn.spatial_shape)),
            np.full((2, *sn.spatial_shape), 0.5),
            _random_composite(sn, np.random.default_rng(1)).boundary,
        )


def test_tail_mismatch_refuses_loudly():
    """[Pattern-4 pin] a cotangent whose spatial-moment tail does not match
    the scheme's raises the typed backstop instead of broadcasting silently
    through the batch VJP (the multi-D mirror of the C2 1-D backstop)."""
    sn = _cart2d_probe_mesh()
    rng = np.random.default_rng(20260802)
    phi = _random_composite(sn, rng)
    # Graft a bogus trailing moment axis onto the DD (tail-less) cotangent
    # (width 2^d = 4, the LD-2D layout the DD scheme does not carry).
    bad_interior = AngularFlux.from_mesh(
        np.asarray(phi.interior.values)[..., None].repeat(4, axis=-1), sn,
        spatial_moments=2,
    )
    bad = FullField(interior=bad_interior, boundary=phi.boundary)
    with pytest.raises(ValueError, match="spatial-moment tail"):
        FullFieldWavefront(sn).loss_action_transpose(_het_sigma(sn, rng), bad)

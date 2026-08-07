r"""The algebra of record for SN sweep acyclicity — when is ``L + C`` triangular?

This module is the **source of truth** for the acyclicity claims in
``docs/theory/foundations/path_integral.rst`` (§"the sweep is a
forward substitution") and ``docs/theory/foundations/boundary_conditions.rst``
(the boundary-cycle discussion).  If a triangularity/cycle claim in the RST
cannot be produced by this module, it must be added here first.

Why the module exists
=====================

The sweep's cheapness rests on a **theorem with hypotheses**: for a fixed
ordinate the cell-dependency digraph is acyclic, so the discrete operator is
triangular in walk order and inverts by one forward substitution.  ORPHEUS
states that position honestly —

    Acyclicity is a property of the **(mesh, closure, boundary) triple**, and
    ORPHEUS certifies it *per case* rather than assuming it.

— but until now there was no **executable** object behind the claim.  The
production sweep order is closed-form index arithmetic
(``level_of = local.sum(axis=0)`` in
:mod:`orpheus.sn.loss_representation.sweep_graph`), deliberately *not* a
topological sort, so nothing in the solver ever builds a digraph and nothing
can detect a cycle.  The invariant the design named for that job —
``assert_cycles_are_declared`` — was never implemented, and the per-law
``BoundaryTraceLaw.creates_sweep_cycle`` ClassVar meant to feed it was read by
no production code and has since been **retired**.  That flag could not have
worked even in principle: whether a face's trace back-edge closes a *cycle*
depends on the whole face configuration (``reflective|vacuum`` is acyclic,
``reflective|reflective`` is not), which a boolean on the boundary *kind*
cannot express.  This module supplies the missing object.

Two digraphs, two theorems — do not conflate them
=================================================

**(1) The CELL digraph, per ordinate — the bulk.**  For a fixed ordinate on a
structured Cartesian mesh with an upwind closure, acyclicity is a *lattice*
theorem:

.. math::

    \Phi(i) \;=\; \sum_a \operatorname{sign}(\Omega_a)\, i_a

is a strict potential that **every** dependency edge increases, so no cycle can
close.  Note what the proof uses — the mesh's lattice order, not the
characteristics.  It is a mesh theorem, exactly as strong as its hypotheses,
and it fails for an unstructured mesh (the Pautz/Plimpton cycle-breaking
problem).  This module does not re-derive it; the production certificate is the
assembled-matrix gate ``triu(PᵀMP, 1) == 0`` in
``tests/sn/sweep/test_assembly_mode.py``.

**(2) The TRACE digraph, across ordinates — the boundary.**  This is the one
that decides whether the *whole* within-group problem sweeps in one pass, and
it is the one this module builds.  Its nodes are trace slots
:math:`(\text{face}, \text{ordinate})`; its edges are

* a **sweep edge** ``inflow(n) → outflow(n)`` per ordinate — the transport
  carrying information across the domain; and
* a **boundary edge** per reflecting/periodic face — the law feeding an
  outflow slot back into an inflow slot.

The folklore this refutes
=========================

    "A reflective boundary creates a cycle, so the reflected coupling must be
    extracted from the sweep."

**False as stated.**  A *single* reflecting face is a forward edge: the sweep
visits :math:`\mu < 0` before :math:`\mu > 0`, so the mirror feeds information
strictly downstream.  ORPHEUS relies on this — the curvilinear :math:`r = 0`
pole mirror is a specular reflection kept **inside** the walk.  A cycle needs a
closed **loop**, e.g. both faces of a slab reflecting.  The honest criterion is
therefore not a boolean on the boundary *kind* but a strongly-connected-
component decomposition of the trace digraph — which is what
:func:`derive_slab_trace_acyclicity` computes.

The consequence for naming, recorded because it is easy to lose: triangularity
is a property of a *configuration*, never of the operator alone.  An operator
name may state what the object **is** (streaming + collision); it must not
assert sweepability, which a new mesh or boundary can falsify.

References
----------
* Adams & Larsen (2002), *Fast iterative methods for discrete-ordinates
  particle transport calculations*, §VI (parallel sweeps).
* Pautz (2002), *An Algorithm for Parallel S_n Sweeps on Unstructured Meshes*
  — the cycle-breaking problem this module's criterion generalises to.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

__all__ = [
    "TraceDigraph",
    "TraceSlot",
    "build_slab_trace_digraph",
    "derive_slab_trace_acyclicity",
]

#: Boundary kinds this module models, by their trace-coupling structure.
#: ``vacuum`` / ``prescribed_inflow`` supply the inflow as DATA (no edge);
#: ``reflective`` maps an outflow to the mirror ordinate's inflow on the SAME
#: face; ``white`` couples EVERY outflow on the face to EVERY inflow on it;
#: ``periodic`` maps an outflow to the SAME ordinate's inflow on the OPPOSITE
#: face.
BoundaryKind = Literal["vacuum", "prescribed_inflow", "reflective", "white", "periodic"]

_FACES: tuple[str, str] = ("left", "right")


@dataclass(frozen=True, order=True)
class TraceSlot:
    """One node of the trace digraph: an ordinate's flux on one face.

    ``role`` is derived, never chosen: for a slab the face an ordinate ENTERS
    through is fixed by :math:`\\operatorname{sign}(\\mu)`.
    """

    face: str
    ordinate: int
    role: Literal["inflow", "outflow"]

    def __str__(self) -> str:  # pragma: no cover - display only
        return f"{self.role}({self.face}, n={self.ordinate})"


@dataclass(frozen=True)
class TraceDigraph:
    """The (face, ordinate) trace-dependency digraph of a 1-D slab.

    Attributes
    ----------
    slots :
        Node list; ``adjacency[i, j] == 1`` means slot ``j`` depends on slot
        ``i`` (information flows ``i → j``).
    adjacency :
        Dense 0/1 adjacency, ``(n_slots, n_slots)``.
    sweep_edges, boundary_edges :
        The two edge families, kept separate so a cycle can be attributed to
        the boundary law that closed it.
    """

    slots: tuple[TraceSlot, ...]
    adjacency: np.ndarray
    sweep_edges: tuple[tuple[int, int], ...]
    boundary_edges: tuple[tuple[int, int], ...]

    # ── structure ────────────────────────────────────────────────────────

    def strongly_connected_components(self) -> tuple[tuple[int, ...], ...]:
        """SCCs of the digraph, as index tuples, largest first.

        A component of size > 1 is a genuine cycle: no ordering of the slots
        can place all its edges below the diagonal.
        """
        n, labels = connected_components(
            csr_matrix(self.adjacency), directed=True, connection="strong",
        )
        comps = [tuple(np.flatnonzero(labels == k).tolist()) for k in range(n)]
        return tuple(sorted(comps, key=len, reverse=True))

    @property
    def is_acyclic(self) -> bool:
        """``True`` iff every SCC is a single slot — i.e. the trace digraph
        is a DAG and the whole within-group problem sweeps in ONE pass."""
        return all(len(c) == 1 for c in self.strongly_connected_components())

    def cyclic_components(self) -> tuple[tuple[int, ...], ...]:
        """The nontrivial (size > 1) SCCs — the cycles, if any."""
        return tuple(
            c for c in self.strongly_connected_components() if len(c) > 1
        )

    def topological_order(self) -> np.ndarray | None:
        """A sweep order making every edge forward, or ``None`` if cyclic.

        Kahn's algorithm.  The returned permutation is exactly the ordering
        under which :meth:`permuted_dependency_matrix` is strictly lower-triangular —
        i.e. the order a one-pass forward substitution would visit the slots.
        """
        a = self.adjacency
        indeg = a.sum(axis=0).astype(int)
        ready = [i for i in range(len(self.slots)) if indeg[i] == 0]
        order: list[int] = []
        while ready:
            # deterministic: smallest index first, so the order is reproducible
            ready.sort()
            i = ready.pop(0)
            order.append(i)
            for j in np.flatnonzero(a[i]):
                indeg[j] -= 1
                if indeg[j] == 0:
                    ready.append(int(j))
        return np.asarray(order, dtype=np.intp) if len(order) == len(self.slots) else None

    def dependency_matrix(self) -> np.ndarray:
        r"""The digraph in **operator convention**: ``M[i, j] != 0`` means the
        equation for slot ``i`` READS slot ``j``.

        This is the transpose of :attr:`adjacency` (which is edge-convention,
        ``i → j``), and it is the convention every triangularity gate in the
        codebase uses — a discrete operator's row is its equation and its
        columns are the unknowns that equation depends on.  Keeping ONE
        convention across the tree is deliberate: the acyclicity test here is
        then literally the production gate's,
        ``np.triu(permuted, k=1) == 0``
        (``tests/sn/sweep/test_assembly_mode.py::test_g2_walk_order_triangularity_is_exact``).
        """
        return self.adjacency.T

    def permuted_dependency_matrix(
        self, order: np.ndarray | None = None,
    ) -> np.ndarray:
        """:meth:`dependency_matrix` permuted to sweep order — **the matrix to
        look at**.

        Strictly LOWER-triangular (``triu(·, 1) == 0``) ⟺ every equation reads
        only already-solved slots ⟺ the configuration sweeps in ONE pass.
        When the digraph is cyclic no such order exists, and the fallback
        :meth:`mu_negative_first_order` exposes the surviving above-diagonal entries:
        those are exactly the **back-edges** a Gauss-Seidel schedule must lag.
        """
        if order is None:
            order = self.topological_order()
            if order is None:
                order = self.mu_negative_first_order()
        return self.dependency_matrix()[np.ix_(order, order)]

    def mu_negative_first_order(self) -> np.ndarray:
        r"""The FIXED order "all :math:`\mu<0` ordinates, then all
        :math:`\mu>0`" (inflow before outflow within each).

        **This is not a universal sweep order, and the difference is the
        point.**  It is triangular for a LEFT-reflecting slab (the
        :math:`\mu<0` sweep reaches the left face first, so the mirror feeds
        strictly downstream) and NOT triangular for a RIGHT-reflecting one,
        even though both are acyclic — a single reflection is one-pass only
        when the sweep DIRECTION matches the reflecting face.  Acyclicity
        (:meth:`strongly_connected_components`) asks whether SOME order
        exists; this asks whether THIS one works.  Use it as the control that
        separates the two questions, and — on a genuinely cyclic
        configuration — to expose the back-edges a Gauss-Seidel schedule must
        lag.
        """
        idx = {s: i for i, s in enumerate(self.slots)}
        negative = [s for s in self.slots if s.face == "right" and s.role == "inflow"]
        order: list[int] = []
        for group in (negative, [s for s in self.slots
                                if s.face == "left" and s.role == "inflow"]):
            for slot in sorted(group, key=lambda s: s.ordinate):
                order.append(idx[slot])
                partner = next(
                    t for t in self.slots
                    if t.ordinate == slot.ordinate and t.role == "outflow"
                )
                order.append(idx[partner])
        return np.asarray(order, dtype=np.intp)

    def describe_cycles(self) -> list[str]:
        """Human-readable attribution: which boundary edges close each cycle."""
        out: list[str] = []
        bset = set(self.boundary_edges)
        for comp in self.cyclic_components():
            members = set(comp)
            closing = [
                f"{self.slots[i]} -> {self.slots[j]}"
                for (i, j) in sorted(bset)
                if i in members and j in members
            ]
            out.append(
                f"SCC of {len(comp)} slots "
                f"[{', '.join(str(self.slots[i]) for i in comp)}]"
                + (f"; closed by boundary edges: {'; '.join(closing)}" if closing else "")
            )
        return out


def build_slab_trace_digraph(
    mu: np.ndarray,
    reflection_index: np.ndarray,
    *,
    bc_left: BoundaryKind,
    bc_right: BoundaryKind,
) -> TraceDigraph:
    r"""Assemble the trace digraph for a 1-D slab.

    Parameters
    ----------
    mu :
        Axis cosines :math:`\mu_n`, shape ``(N,)``.  Only the SIGN is used —
        acyclicity is a combinatorial property, not a numerical one.
    reflection_index :
        ``reflection_index[n]`` is the partner ordinate under reflection across
        the slab axis — the ``.indices`` of
        ``Quadrature.ordinate_permutation(RigidMotion.reflection(normal=ê_x))``.
    bc_left, bc_right :
        The boundary kinds.  See :data:`BoundaryKind`.

    Notes
    -----
    Ordinates with :math:`\mu_n = 0` do not stream along the axis; they carry
    no trace edge (the production walk short-circuits them with
    :math:`Q/\Sigma_t`) and appear as isolated slots.
    """
    mu = np.asarray(mu, dtype=float)
    reflection_index = np.asarray(reflection_index, dtype=int)
    n_ord = mu.size

    def entry_face(n: int) -> str:
        return "left" if mu[n] > 0.0 else "right"

    def exit_face(n: int) -> str:
        return "right" if mu[n] > 0.0 else "left"

    slots: list[TraceSlot] = []
    for n in range(n_ord):
        if mu[n] == 0.0:
            continue
        slots.append(TraceSlot(entry_face(n), n, "inflow"))
        slots.append(TraceSlot(exit_face(n), n, "outflow"))
    slots_t = tuple(slots)
    idx = {s: i for i, s in enumerate(slots_t)}

    a = np.zeros((len(slots_t), len(slots_t)), dtype=int)
    sweep: list[tuple[int, int]] = []
    bnd: list[tuple[int, int]] = []

    # ── sweep edges: inflow --(transport across the slab)--> outflow ──
    for n in range(n_ord):
        if mu[n] == 0.0:
            continue
        i = idx[TraceSlot(entry_face(n), n, "inflow")]
        j = idx[TraceSlot(exit_face(n), n, "outflow")]
        a[i, j] = 1
        sweep.append((i, j))

    # ── boundary edges: outflow --(the law)--> inflow ──
    for face, kind in (("left", bc_left), ("right", bc_right)):
        if kind in ("vacuum", "prescribed_inflow"):
            continue  # inflow is GIVEN data — the digraph has no in-edge
        outs = [s for s in slots_t if s.face == face and s.role == "outflow"]
        ins = [s for s in slots_t if s.face == face and s.role == "inflow"]
        for s in outs:
            if kind == "reflective":
                partner = int(reflection_index[s.ordinate])
                tgt = TraceSlot(face, partner, "inflow")
                if tgt in idx:
                    a[idx[s], idx[tgt]] = 1
                    bnd.append((idx[s], idx[tgt]))
            elif kind == "white":
                # couples EVERY outflow on the face to EVERY inflow on it
                for t in ins:
                    a[idx[s], idx[t]] = 1
                    bnd.append((idx[s], idx[t]))
            elif kind == "periodic":
                other = "right" if face == "left" else "left"
                tgt = TraceSlot(other, s.ordinate, "inflow")
                if tgt in idx:
                    a[idx[s], idx[tgt]] = 1
                    bnd.append((idx[s], idx[tgt]))
            else:  # pragma: no cover - guarded by the Literal
                raise ValueError(f"unmodelled boundary kind {kind!r}")

    return TraceDigraph(slots_t, a, tuple(sweep), tuple(bnd))


def derive_slab_trace_acyclicity(
    mu: np.ndarray, reflection_index: np.ndarray,
) -> dict[tuple[str, str], TraceDigraph]:
    """Build the digraph for the canonical slab BC combinations.

    Returns a mapping ``(bc_left, bc_right) -> TraceDigraph`` for the four
    cases that settle the folklore:

    ============================  =========  ==============================
    configuration                 acyclic?   why
    ============================  =========  ==============================
    vacuum   | vacuum             YES        no boundary edge at all
    reflective | vacuum           **YES**    ONE reflecting face is a
                                             FORWARD edge — the refutation
    vacuum | reflective           **YES**    mirror image of the above
    reflective | reflective       NO         a closed LOOP: the lattice cell
    ============================  =========  ==============================
    """
    cases: tuple[tuple[BoundaryKind, BoundaryKind], ...] = (
        ("vacuum", "vacuum"),
        ("reflective", "vacuum"),
        ("vacuum", "reflective"),
        ("reflective", "reflective"),
    )
    return {
        (bl, br): build_slab_trace_digraph(
            mu, reflection_index, bc_left=bl, bc_right=br,
        )
        for bl, br in cases
    }

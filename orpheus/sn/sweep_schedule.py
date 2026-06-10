r"""Octant-group sweep SCHEDULE — the polymorphic Jacobi / Gauss-Seidel strategy.

Phase 3 (SI Gauss-Seidel rate recovery). The Wave O BC extraction made the
transport sweep BARE and applies the reflective coupling ``B`` externally, which
turned the intra-sweep Gauss-Seidel reflective coupling into inter-sweep Jacobi
(``B`` fully lagged) — same converged fixed point, slower SI spectral rate. The
recovery interleaves the external ``−B`` reflect at octant-group granularity
inside the SI resolvent, realizing the forward substitution ``(L+C−B_lower)⁻¹``.

Jacobi and Gauss-Seidel are the SAME uniform sweep-and-reflect loop differing
ONLY in this schedule (there is NO ``if jacobi/gs`` branch in the iteration — the
splitting is selected ONCE, by choosing the schedule):

* **Jacobi** — ONE group containing every octant, with NO inter-group reflect.
  All octants read the same frozen inflow seed (``B·ψₙ``); identical to the
  pre-recovery bare all-octants sweep.
* **Gauss-Seidel** — one group per in-plane octant (:class:`OctantLabel`), in
  quadrature sweep order; after each group its reflective OUTGOING faces are
  re-reflected (the face-restricted
  :meth:`~orpheus.sn.boundary_operator.SNBoundaryOperator.reflect_into_inflow`),
  so a LATER group reads the fresh current-iterate reflection. Octants swept
  before their specular partner keep the lagged seed (the cyclic ``B_upper``
  back-edges — e.g. a both-faces-reflective axis is a 2-cycle ⟹ only PARTIAL
  one-pass G-S); octants swept after read the fresh value (the order-respecting
  ``B_lower`` edges).

The schedule is a **mesh-time derived object** — it depends only on the
quadrature's octant partition + the mesh's reflective-face set, not on fluxes /
sources / iteration state — so it is built once and reused across every SI
iterate (the same lifetime contract as
:class:`~orpheus.sn.sweep_graph.SweepDependencyGraph`).

The fixed point is INVARIANT under the schedule: any consistent splitting of
``(L+C−S−B)ψ=q`` shares the dominant solution ψ\* (at convergence the seed and
all re-reflects equal ``B·ψ\*``). The schedule changes only the SI spectral rate
(``ρ_J = c`` Jacobi vs ``ρ_GS ≈ c²`` for the symmetric reflective model problem).
Krylov is splitting-invariant and ignores the schedule entirely.

See also
========

* :mod:`orpheus.sn.sweep_graph` — :class:`OctantLabel` + the per-octant
  :class:`SweepDependencyGraph` (the cell DAG the schedule's groups are swept
  through; a distinct concept — cell causality vs octant ordering).
* ``.claude/agent-memory/explorer/si_gs_substep3_carve_substrate.md`` — the
  carve substrate map (octant↔face geometry, the resolvent solve path).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from orpheus.sn.sweep_graph import OctantLabel

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


__all__ = ["OctantSweep", "OctantSweepGroup", "SweepSchedule"]


# An octant streaming in the +axis direction outflows through the axis-max face;
# a grazing (sign 0) axis has NO net outflow on that axis (omitted below).
_OUT_FACE = {
    ("x", +1): "xmax", ("x", -1): "xmin",
    ("y", +1): "ymax", ("y", -1): "ymin",
}


@dataclass(frozen=True, slots=True)
class OctantSweep:
    """One octant's sweep unit: its in-plane :class:`OctantLabel` (selects the
    per-octant :class:`~orpheus.sn.sweep_graph.SweepDependencyGraph`) + the
    ordinate indices into the ``(N, …)`` angular axis it sweeps."""

    label: OctantLabel
    indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class OctantSweepGroup:
    """Octants swept together over ONE frozen inflow; after the group, its
    ``reflect_faces`` are re-reflected so later groups read the fresh inflow.

    ``reflect_faces`` is empty for the Jacobi group, for a grazing / pure-z
    octant (no in-plane outflow), and for non-reflective (vacuum / white)
    boundaries — in all of which the inter-group reflect is a no-op.
    """

    sweeps: tuple[OctantSweep, ...]
    reflect_faces: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SweepSchedule:
    """Polymorphic octant-group schedule for the SI resolvent's uniform
    sweep-and-reflect loop. See the module docstring for the Jacobi/G-S split.
    """

    groups: tuple[OctantSweepGroup, ...]
    kind: str  # "jacobi" | "gauss_seidel" — diagnostic / introspection only

    @classmethod
    def jacobi(cls, sn_mesh: "SNMesh") -> "SweepSchedule":
        """One group, all octants, no inter-group reflect — the bare all-octants
        sweep with the whole ``B·ψₙ`` seed frozen for the entire sweep."""
        sweeps = tuple(_octant_sweep(entry) for entry in sn_mesh.quad.octants)
        return cls(
            groups=(OctantSweepGroup(sweeps=sweeps, reflect_faces=()),),
            kind="jacobi",
        )

    @classmethod
    def gauss_seidel(cls, sn_mesh: "SNMesh") -> "SweepSchedule":
        """One group per in-plane octant, in quadrature sweep order; each group
        re-reflects the reflective faces its octants outflow through.

        Octant partition entries that share an in-plane :class:`OctantLabel`
        (they differ only in the out-of-plane sign_z the 2-D sweep ignores) are
        MERGED into one group so a face's outflow is complete before it is
        reflected.
        """
        reflective = _reflective_faces(sn_mesh)
        ordered: list[OctantLabel] = []
        by_label: dict[OctantLabel, list[OctantSweep]] = {}
        for entry in sn_mesh.quad.octants:
            sweep = _octant_sweep(entry)
            if sweep.label not in by_label:
                by_label[sweep.label] = []
                ordered.append(sweep.label)
            by_label[sweep.label].append(sweep)

        # Assign each reflective face to the LAST in-plane octant group (in
        # sweep order) that OUTFLOWS through it — reflecting only after that
        # group guarantees the face's outflow is COMPLETE (every octant that
        # streams out through it has been swept this pass), so the reflected
        # inflow is consistent.
        #
        # ⚠ Correctness (NOT just rate): a face is shared by EVERY octant whose
        # sign on that axis matches (e.g. ``xmax`` ← all +x octants).  For an
        # axis-aligned quadrature (``product`` — single-face octants) each face
        # has exactly ONE outflowing group, so "last" = that group.  But for a
        # diagonal / spherical cubature (``lebedev`` / ``level_symmetric`` —
        # each octant outflows TWO faces) a face is shared by ≥2 groups;
        # reflecting after the FIRST would absorb the not-yet-swept octants'
        # SEED value (the wavefront is rebuilt + seeded each solve, so their
        # outflow slots still hold the inflow seed, NOT real outflow) and
        # reflect garbage — converging to a WRONG fixed point.  Deferring to
        # the LAST outflowing group fixes this; octants reading the face that
        # are swept BEFORE its reflect keep the lagged seed (the cyclic
        # back-edge → partial one-pass G-S), which is always valid.
        last_group_for_face: dict[str, int] = {}
        for gi, label in enumerate(ordered):
            for f in _outgoing_faces(label):
                if f in reflective:
                    last_group_for_face[f] = gi   # later gi wins → the last
        reflect_by_group: dict[int, list[str]] = {
            gi: [] for gi in range(len(ordered))
        }
        for face, gi in last_group_for_face.items():
            reflect_by_group[gi].append(face)

        groups = tuple(
            OctantSweepGroup(
                sweeps=tuple(by_label[label]),
                reflect_faces=tuple(sorted(reflect_by_group[gi])),
            )
            for gi, label in enumerate(ordered)
        )
        return cls(groups=groups, kind="gauss_seidel")


def _octant_sweep(entry) -> OctantSweep:
    """Project a quadrature octant partition entry to its in-plane sweep unit.

    The partition ``.label`` is ``(sign_x[, sign_y[, sign_z]])`` (∈ {−1, 0, +1});
    the 2-D Cartesian sweep is invariant under sign_z, so the in-plane
    :class:`OctantLabel` drops it (mirroring ``_sweep_2d_wavefront``).
    """
    label = entry.label
    sign_x = int(label[0])
    sign_y = int(label[1]) if len(label) >= 2 else 0
    return OctantSweep(
        label=OctantLabel((sign_x, sign_y)),
        indices=tuple(int(i) for i in entry.indices),
    )


def _outgoing_faces(label: OctantLabel) -> tuple[str, ...]:
    """The boundary faces an octant OUTFLOWS through (strict sign — a grazing
    ``sign == 0`` axis contributes no net outflow on that axis)."""
    faces: list[str] = []
    for axis, sign in (("x", label.sign_x), ("y", label.sign_y)):
        face = _OUT_FACE.get((axis, sign))
        if face is not None:
            faces.append(face)
    return tuple(faces)


def _reflective_faces(sn_mesh: "SNMesh") -> frozenset[str]:
    """The mesh's SPECULAR-reflective boundary faces.

    Vacuum (no coupling) and white (couples ALL ordinates on the face ⟹ the
    octant-order G-S degenerates to Jacobi) are EXCLUDED — only specular
    reflection admits the order-respecting forward-substitution acceleration.
    """
    return frozenset(
        face
        for face in sn_mesh.trace.layout.faces
        if getattr(sn_mesh, f"bc_{face}") == "reflective"
    )

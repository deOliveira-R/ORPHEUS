r"""The boundary-trace function space :math:`\Gamma = \partial\Omega \times S^2`.

A *trace space* is the :class:`FunctionSpace` that lives on the
boundary :math:`\partial\Omega` of the spatial domain, carrying the
angular degrees of freedom on each boundary face. It is the
domain/codomain space of every boundary operator (vacuum, albedo,
reflective, white, prescribed-inflow) and the storage space of the
boundary flux.

One space, two directional *selectors* — not three types
========================================================

Issues #205 / #201 (the View-G field-vocabulary refactor) collapse the
three previously-separate boundary-space notions into **one** concrete
:class:`TraceSpace`:

* the per-face ``InflowTraceSpace`` / ``OutflowTraceSpace`` pair
  (Wave 2, ``transient-giggling-cake``), and
* the ad-hoc ``FunctionSpace("sn_boundary_flat")`` that
  :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` built
  for its flat storage, and
* the dead ``boundary_trace_space()`` factory.

The unification rests on one observation: **inflow and outflow are
operations on a single space, not two spaces.** Whether an ordinate is
incoming or outgoing at a face is a *predicate* — :math:`\mathrm{sign}
(\Omega \cdot \hat n_f)` — evaluated against the same trace data, not a
property of the space's identity. So :class:`TraceSpace` stores the
*signed* projection :math:`\Omega \cdot \hat n_f` once, per face, and
exposes :meth:`inflow_indices_for_face` / :meth:`outflow_indices_for_face`
as selectors over it. (#208 will promote these to projection
*operators*; the :math:`|\Omega\cdot\hat n|`-weighted boundary inner
product they live in is now installed — see below.)

Whole-boundary storage + per-face access
=========================================

The space is the **whole** boundary: ``shape == (layout.total_size,)``
where ``layout`` is the :class:`~orpheus.numerics.face_layout.FaceLayout`
that packs every face into one flat buffer (the same descriptor
:class:`BoundaryFlux` consumes). Per-face access is via the layout's
slot (``layout.faces[face].slice_view``) plus the per-face row of the
signed-projection table; the old per-face ``(N, ng)`` "space" is now a
*derived view*, not a class.

Inner product is the **partial-current metric**
:math:`G_s = |\Omega\cdot\hat n_f|\odot w_n` (the cosine-weighted angular
quadrature), installed at construction by :func:`_build_trace_metric_weights`.
This is the physically-correct boundary inner product under which the
``BoundaryOperator`` Hilbert adjoints (``B.H``) — reflective and white —
are correct (Wave O / O.2b, #208). The metric is group-independent (a
weight in angle, not energy). It is read ONLY by the adjoint path
(:class:`~orpheus.numerics.operator._AdjointOperator` and
:meth:`FunctionSpace.inner_product`); the forward sweep/matvec never reads
it, so installing it leaves every forward result bit-identical. (Before
O.2b this slot was Euclidean ``None``, matching the legacy
``sn_boundary_flat`` storage space.)

Geometric convention
====================

For each face with outward unit normal :math:`\hat n_f`, the signed
projection :math:`\Omega_n \cdot \hat n_f = \mathrm{sign}_f \cdot
\mu_{\text{axis}(f)}` classifies ordinate :math:`\Omega_n` as:

* **Inflow** iff :math:`\Omega_n \cdot \hat n_f < -\epsilon`
  (direction points INTO the domain),
* **Outflow** iff :math:`\Omega_n \cdot \hat n_f > +\epsilon`
  (direction points OUT of the domain),
* **Tangential** iff :math:`|\Omega_n \cdot \hat n_f| \leq \epsilon`
  (grazes the face; in NEITHER selector).

Principled tolerance (not a magic number)
-----------------------------------------

``_TANGENTIAL_EPS = 4 * np.finfo(np.float64).eps`` (:math:`\approx
8.9\times10^{-16}`). It is a safety factor over the IEEE-754
dot-product round-off bound :math:`d \cdot u` (:math:`d \leq 3` spatial
dimensions, :math:`u = \epsilon_{\mathrm{mach}}/2` the unit round-off)
for the unit-vector projection :math:`\Omega\cdot\hat n = \langle \hat
n, \mu\rangle`. Empirically (``eps_probe.py``, Gauss-Legendre
``N=2..64`` + Lebedev orders ``3..53``): nominally-tangential cosines
are **exactly** ``0.0`` (quadrature symmetry — odd-N central node, all
off-axis 1-D components, Lebedev axis nodes), while the smallest
*genuine* cosine is :math:`2.44\times10^{-2}`. The gap
:math:`[0, 0.024]` spans ~14 orders, so this eps sits 4× above the
round-off floor and :math:`2.7\times10^{13}\times` below any genuine
projection — making the inflow/outflow masks **bit-identical** to both
the operator's former ``1e-15`` and the realizer's former ``1e-12``
(the band ``(eps, 1e-12)`` is empty). :func:`test_eps_below_min_genuine_cosine`
guards the gap so a future quadrature cannot silently violate it.

Coord-system coverage
=====================

The face → outward-normal table is keyed on the **layout's** face
names (``SNMesh.boundary_face_layout``), which is the single source of
truth for "which faces exist":

* **Mesh1D slab** — two faces ``xmin`` / ``xmax``; outward normals
  :math:`\mp\hat x`.
* **Mesh1D curvilinear** (sphere / cylinder) — ONE face ``xmax`` (the
  outer radius :math:`r=R`). The geometric pole :math:`r=0` is a
  *regularity/symmetry condition* handled by the angular sweep's
  pole closure (``MorelMontryAngularSweep`` / Carlson seed), **not** a
  boundary face — there is no surface there and no inflow to impose.
  This is why the curvilinear layout has no inner face (it is NOT a
  ``left/right`` pair): a solid sphere has exactly one boundary.
* **Mesh2D Cartesian** — four faces ``xmin`` / ``xmax`` / ``ymin`` /
  ``ymax``; a 3-axis Cartesian mesh (C5.5, #225) all six. 2-D
  cylindrical ``(r, z)`` has no SN sweep and cannot become an SNMesh
  (refused at the axis conversion during construction — since C5.3 the
  trace itself is geometry-blind and never sees a mesh).

References
----------

* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
  Transport*. American Nuclear Society. §3.7 (boundary trace operators
  in the discrete-ordinates setting), §6 (curvilinear angular
  redistribution / starting-direction closure at :math:`r=0`).
* ``.claude/plans/field_role_typing_view_g.md`` — A.2/A.3 TraceSpace
  unification design (View-G, signed-:math:`\Omega\cdot\hat n`,
  principled eps, face-naming reconciliation).
* Issue #208 comment (2026-05-31) — the
  :math:`|\Omega\cdot\hat n|`-weighted partial-current boundary inner
  product for the Wave-O adjoint work (landed Phase 4 / O.2b).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.face_layout import AXIS_NAMES
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.numerics.quadrature import Quadrature


__all__ = ["TraceSpace"]


# Tangential tolerance for the unit-vector projection ``Ω · n``. A
# safety factor (×4) over the IEEE-754 dot-product round-off bound for a
# 3-component unit-vector inner product. Empirically bit-identical to
# the legacy ``1e-15`` (operator) and ``1e-12`` (realizer) tolerances;
# see the module docstring + :func:`test_eps_below_min_genuine_cosine`.
_TANGENTIAL_EPS: float = 4.0 * np.finfo(np.float64).eps


# ─────────────────────────────────────────────────────────────────────
# Face → outward-normal table
# ─────────────────────────────────────────────────────────────────────
#
# Each entry maps a face name to ``(axis_index, sign)`` where
# ``axis_index`` selects ``(mu_x, mu_y, mu_z)[axis_index]`` and ``sign``
# is the outward-normal sign (±1). The signed projection is
# ``Ω · n = sign * mu[axis]``. This single table serves every supported
# mesh: 1-D meshes use ``xmin`` / ``xmax`` (radial axis is the x-axis;
# curvilinear has ``xmax`` only), 2-D Cartesian the four x/y faces,
# 3-D Cartesian all six.
#
# C5.3 (#225): DERIVED from :data:`~orpheus.numerics.face_layout.AXIS_NAMES`
# — the same crosswalk every face-name producer renders through
# (``FaceLabel.face_name``, the sweep schedule, the walk) — instead of a
# hand-listed 4-face transcription that silently lacked the z faces.
# ``min`` is the −axis face (outward normal −ê_axis), ``max`` the +axis.

_FACE_NORMALS: dict[str, tuple[int, int]] = {
    f"{name}{suffix}": (axis, sign)
    for axis, name in enumerate(AXIS_NAMES)
    for suffix, sign in (("min", -1), ("max", +1))
}


def _quadrature_axis(quadrature: "Quadrature", axis: int) -> np.ndarray:
    """Return ``mu_x`` for axis=0, ``mu_y`` for axis=1, ``mu_z`` for axis=2.

    Falls back to a zero array when the requested axis is not present on
    the quadrature object — this lets a 1-D :class:`Quadrature` (which
    has ``mu_x`` populated and ``mu_y == 0``) feed the same predicate
    logic without special-casing.
    """
    mu_x = np.asarray(quadrature.mu_x)
    if axis == 0:
        return mu_x
    if axis == 1:
        return np.asarray(getattr(quadrature, "mu_y", np.zeros_like(mu_x)))
    if axis == 2:
        return np.asarray(getattr(quadrature, "mu_z", np.zeros_like(mu_x)))
    raise ValueError(f"axis must be 0, 1, or 2; got {axis}")


def _build_omega_dot_n(
    quadrature: "Quadrature",
    faces: tuple[str, ...],
) -> NDArray:
    r"""Build the signed projection table :math:`\Omega \cdot \hat n_f`.

    Returns a ``(n_faces, n_ordinates)`` float array whose row ``f`` is
    :math:`\mathrm{sign}_f \cdot \mu_{\text{axis}(f)}` — the outward
    projection of every ordinate onto face ``f``'s normal. Inflow /
    outflow / tangential are derived from its sign on demand.

    C5.3 (#225): geometry-blind — every datum comes from the quadrature
    and the face NAMES (the axis-aligned outward normals are implied by
    the ``"{axis}{min|max}"`` convention). The former ``mesh``
    parameter was gate-only: its curvilinear-``Mesh2D`` refusal is
    unreachable (such a mesh cannot become an ``SNMesh`` — the axis
    conversion at construction refuses it), and the isinstance check
    carried no data.
    """
    n_ord = int(quadrature.N)
    omega_dot_n = np.zeros((len(faces), n_ord), dtype=float)
    for f_idx, face_name in enumerate(faces):
        try:
            axis, sign = _FACE_NORMALS[face_name]
        except KeyError as exc:
            raise ValueError(
                f"Unknown face name {face_name!r}; valid faces are "
                f"{sorted(_FACE_NORMALS)}"
            ) from exc
        # C5.5 (#225) fail-loud: a layout naming an axis-k FACE demands
        # GENUINE mu_k on the quadrature. Discriminate on VALUE, not
        # attribute presence — the per-axis cosines are properties that
        # zero-pad past the cubature's intrinsic dimensionality (e.g.
        # 1-D Gauss-Legendre carries mu_z == zeros(N), never an absent
        # attribute), so an attribute test can never fire. A boundary
        # face whose normal-axis cosines are ALL zero has Ω·n ≡ 0 for
        # every ordinate — a rank-mismatch (a z face on a quadrature
        # with no third cosine) that zero-padding would silently
        # misclassify as all-tangential (neither inflow nor outflow).
        mu_axis = _quadrature_axis(quadrature, axis)
        if not np.any(mu_axis):
            raise ValueError(
                f"Face {face_name!r} requires genuine "
                f"mu_{AXIS_NAMES[axis]} cosines, but every ordinate of "
                f"the quadrature has mu_{AXIS_NAMES[axis]} == 0 — a "
                f"rank-mismatch between the face layout and the "
                f"angular cubature."
            )
        omega_dot_n[f_idx] = sign * mu_axis
    return omega_dot_n


def _build_trace_metric_weights(
    omega_dot_n: NDArray,
    quad_weights: NDArray,
    layout: "FaceLayout",
) -> NDArray:
    r"""Build the partial-current boundary metric :math:`G_s = |\Omega\cdot\hat n_f|\odot w_n`.

    The trace Hilbert metric is the **partial current** weight: pairing two
    boundary fields contracts angle against :math:`|\Omega\cdot\hat n_f|\,w_n`,
    i.e. the cosine-weighted angular quadrature (Lewis & Miller §3.7; the
    boundary inner product under which reflective/white BCs are self-adjoint).
    Wave O / O.2b (#208) — replaces the legacy Euclidean (``None``) metric.

    Returns the flat ``(layout.total_size,)`` diagonal-weight array that
    :meth:`FunctionSpace.inner_product` broadcasts against the trace state.
    The metric is **purely angular** — :math:`|\Omega\cdot\hat n_f|\,w_n`
    depends only on the ordinate (axis 0 of every face slot), not on energy
    group or on spatial position along the face. So for a face slot of shape
    ``(N, ng)`` (1-D) or ``(N, ng, n_face_cells)`` (2-D edge) the ``(N,)``
    cosine weight is broadcast across ALL trailing (group / spatial) axes.

    The row order of ``omega_dot_n`` matches ``tuple(layout.faces)`` (both
    derive from the same ordered layout), so ``enumerate(layout.faces)``
    aligns face slots with projection rows.
    """
    weights_flat = np.zeros((int(layout.total_size),), dtype=float)
    w_n = np.asarray(quad_weights, dtype=float)  # (N,)
    for f_idx, face_name in enumerate(layout.faces):
        slot = layout.faces[face_name]
        face_w = np.abs(omega_dot_n[f_idx]) * w_n  # (N,) = |Ω·n_f| · w_n
        # Ordinate is axis 0 of the slot; reshape to (N, 1, 1, …) so the
        # per-ordinate cosine weight broadcasts across every trailing axis
        # (group, and — in 2-D — the cells along the boundary edge).
        face_w_axis0 = face_w.reshape((face_w.shape[0],) + (1,) * (len(slot.shape) - 1))
        flat_face = np.broadcast_to(face_w_axis0, slot.shape).reshape(-1)
        weights_flat[slot.offset : slot.offset + slot.flat_size] = flat_face
    return weights_flat


# ─────────────────────────────────────────────────────────────────────
# The trace space
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class TraceSpace(FunctionSpace):
    r"""The boundary-trace function space (View-G, role-agnostic).

    One concrete space for the whole boundary :math:`\Gamma`. Inflow and
    outflow are *selectors* over the signed projection
    :math:`\Omega\cdot\hat n`, not separate types. See the module
    docstring for the unification rationale and geometric convention.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`FunctionSpace`. ``name`` is ``"sn_trace"``
        and ``shape`` is the whole-boundary flat shape
        ``(layout.total_size,)``. ``inner_product_weights`` is the
        partial-current metric :math:`G_s = |\Omega\cdot\hat n_f|\odot w_n`
        (built by :func:`_build_trace_metric_weights`; see the module
        docstring) — NOT Euclidean.
    layout : FaceLayout
        The flat-buffer descriptor (which faces exist, per-face shapes,
        offsets). Carried as ``compare=False`` leaf-data so it does not
        pollute the ``(name, shape)`` identity — two trace spaces on
        meshes of the same total boundary size compare equal regardless
        of their face decomposition. Ordered iteration of
        ``layout.faces`` defines the row order of :attr:`omega_dot_n`.
    omega_dot_n : NDArray, shape ``(n_faces, n_ordinates)``
        The signed projection :math:`\Omega\cdot\hat n_f` per face. The
        single source of truth the inflow/outflow selectors AND the
        operator-side directional masks both read.
    """

    layout: Optional["FaceLayout"] = field(
        default=None, repr=False, compare=False,
    )
    omega_dot_n: Optional[NDArray] = field(
        default=None, repr=False, compare=False,
    )

    @property
    def face_names(self) -> tuple[str, ...]:
        """Ordered face names (matching :attr:`omega_dot_n` row order)."""
        if self.layout is None:
            return ()
        return tuple(self.layout.faces)

    @classmethod
    def from_quadrature_and_layout(
        cls,
        quadrature: "Quadrature",
        layout: "FaceLayout",
    ) -> "TraceSpace":
        r"""Build the trace space from a quadrature and a face layout.

        C5.3 (#225): geometry-blind — the former ``mesh`` parameter
        (``from_mesh_and_quadrature``) was gate-only and is retired;
        every datum comes from the quadrature and the layout's face
        names (axis-aligned outward normals implied by the
        ``"{axis}{min|max}"`` convention).

        Parameters
        ----------
        quadrature : Quadrature
            Angular quadrature exposing ``mu_x`` (always) and ``mu_y`` /
            ``mu_z`` (when applicable).
        layout : FaceLayout
            The boundary face layout (canonically
            :attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout` —
            the single source of truth for which faces exist and their
            flat packing). Its ordered faces drive the
            :attr:`omega_dot_n` rows; its ``total_size`` sets the space
            shape.

        Raises
        ------
        ValueError
            If ``layout`` names a face absent from the normal table.
        """
        faces = tuple(layout.faces)
        omega_dot_n = _build_omega_dot_n(quadrature, faces)
        # Wave O / O.2b (#208): the partial-current boundary metric
        # G_s = |Ω·n_f| ⊙ w_n — the cosine-weighted angular quadrature under
        # which the BoundaryOperator Hilbert adjoints (B.H) are physically
        # correct.  Group-independent; built once at the producer (Pattern 7).
        inner_product_weights = _build_trace_metric_weights(
            omega_dot_n, quadrature.weights, layout,
        )
        return cls(
            name="sn_trace",
            shape=(int(layout.total_size),),
            inner_product_weights=inner_product_weights,
            layout=layout,
            omega_dot_n=omega_dot_n,
        )

    # ── Directional selectors ────────────────────────────────────────

    def _face_row(self, face: str) -> int:
        """Return the :attr:`omega_dot_n` row index for ``face``."""
        if self.omega_dot_n is None or self.layout is None:
            raise RuntimeError(
                "TraceSpace was constructed without omega_dot_n / layout; "
                "use TraceSpace.from_quadrature_and_layout() instead of "
                "the bare dataclass constructor."
            )
        try:
            return self.face_names.index(face)
        except ValueError as exc:
            raise ValueError(
                f"Unknown face {face!r}; available: {self.face_names}"
            ) from exc

    def inflow_indices_for_face(self, face: str) -> np.ndarray:
        r"""Ordinate indices that are inflow at ``face``.

        Inflow iff :math:`\Omega\cdot\hat n_f < -\epsilon` (direction
        points into the domain). Tangential ordinates are excluded.
        """
        row = self.omega_dot_n[self._face_row(face)]  # type: ignore[index]
        return np.flatnonzero(row < -_TANGENTIAL_EPS)

    def outflow_indices_for_face(self, face: str) -> np.ndarray:
        r"""Ordinate indices that are outflow at ``face``.

        Outflow iff :math:`\Omega\cdot\hat n_f > +\epsilon` (direction
        points out of the domain). Tangential ordinates are excluded.
        """
        row = self.omega_dot_n[self._face_row(face)]  # type: ignore[index]
        return np.flatnonzero(row > +_TANGENTIAL_EPS)

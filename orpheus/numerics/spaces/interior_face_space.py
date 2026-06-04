r"""The interior face-flux cochain space :math:`C^1_{\rm int}`.

An *interior face space* is the :class:`FunctionSpace` of the SN wavefront
sweep's **interior** cell-face angular fluxes — the per-ordinate values on
every face *inside* the domain that the wavefront propagates across as it
advances from the inflow trace to the outflow trace. It is the interior
sibling of the boundary
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace`: where the trace
space is the boundary 1-cochain :math:`C^1_\partial`, this is the interior
1-cochain :math:`C^1_{\rm int}`, and together they biproduct-decompose the
full face cochain :math:`C^1 = C^1_{\rm int} \oplus C^1_\partial` (the
#208 ``V_bulk ⊕ V_boundary`` shape, one locus down at the face level).

Native frame — discrete exterior calculus / cochains
=====================================================

The per-ordinate angular flux on faces is a primal **1-cochain**
:math:`\psi^{(1)}_\Omega \in C^1` (a value per oriented face). The
cell-average :math:`\bar\psi` is a **0-cochain**; the diamond-difference
closure :math:`\psi^{\rm out} = 2\bar\psi - \psi^{\rm in}` is the averaging
map :math:`C^1 \to C^0`. The boundary trace is the pullback :math:`\iota^*`
under :math:`\partial\Omega \hookrightarrow \Omega`; "absorption = identity"
is :math:`\iota^* \circ \iota_* = \mathrm{id}` on the boundary chain. (See
the cross-domain-attacker frame memo
``field_role_typing_faceflux_frames.md`` for the full validation.)

Layout-on-space (A.5), minus the directional selectors
======================================================

Like :class:`TraceSpace`, this space carries its
:class:`~orpheus.numerics.face_layout.FaceLayout` as ``compare=False``
leaf-data (the A.5 "layout-on-space" re-home), so the
:class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux` leaf reads
``self.space.layout`` rather than carrying a separate field attribute.
UNLIKE :class:`TraceSpace`, it carries **no** ``omega_dot_n``: the interior
cochain is the transient off-diagonal of a per-octant triangular factor and
is **flux-only** (single role), so it has no inflow/outflow partition — the
role grid (flux / source / residual) is a 0-cochain (cell) concept that only
the boundary 1-chain inherits, via its BC consistency residual. Inner
product is **Euclidean** (``inner_product_weights=None``) — the interior
cochain is sweep-transient working state, not a Krylov/adjoint space.

Axis-parametric (3-D-ready by construction)
===========================================

The interior layout is built per-axis: ONE face-normal field per active
axis, shape ``(N, ng, *cells_with_that_axis_incremented)`` — 1-D = ``x``
only ``(N, ng, nx+1)``; 2-D = ``x`` ``(N, ng, nx+1, ny)`` + ``y``
``(N, ng, nx, ny+1)``; 3-D adds ``z``. :meth:`interior_layout` takes the
axis count as a parameter, so 3-D is a parameter (``axes=(0,1,2)``), not a
new code path (plan ``wavefront_flux_foundation.md`` §3a′ — the foundation
for the future ``nd_foundation.md`` session). No ``Mesh3D`` is reachable
from :meth:`from_mesh` yet (it dispatches 1-D vs 2-D on ``mesh.ny``); the
TYPE is forward-built, the implementation stays 1-D + 2-D.

Unify-after-two
===============

:class:`TraceSpace` (boundary, instance 1) and :class:`InteriorFaceSpace`
(interior, instance 2) share the flat-buffer + :class:`FaceLayout` +
``face_names`` machinery. Per ``feedback_unify_after_two_instances`` the
shared ``FaceSpace`` ABC is lifted only AFTER both concrete instances
exist (a later elegance-enforcer pass) — this is that second instance,
built standalone first.

References
----------

* ``.claude/plans/wavefront_flux_foundation.md`` (§0 cochain frame, §3
  the type + biproduct, §3a′ axis-parametricity).
* ``.claude/agent-memory/cross-domain-attacker/field_role_typing_faceflux_frames.md``
  (native-frame validation: ``C^1 = C^1_int ⊕ C^1_∂``; field+views NOT
  per-face objects; the interior cochain is the off-diagonal of ``L_oct``).
* ``coding-elegance`` Pattern 2 (single source of truth — the layout
  descriptor), Pattern 5 (build the primitive: ``FaceLayout``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


__all__ = ["InteriorFaceSpace"]


#: Per-axis face-normal field names, indexed by axis. The interior cochain
#: has ONE such field per active axis (x-normal, y-normal, z-normal faces).
_AXIS_NAMES: tuple[str, ...] = ("x", "y", "z")


@dataclass(frozen=True)
class InteriorFaceSpace(FunctionSpace):
    r"""The interior face-flux cochain space :math:`C^1_{\rm int}` (View-G).

    One concrete space for ALL interior faces of the SN wavefront sweep, on a
    single flat backing buffer described by an axis-parametric
    :class:`~orpheus.numerics.face_layout.FaceLayout`. See the module
    docstring for the cochain frame and the layout-on-space rationale.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited from :class:`FunctionSpace`. ``name`` is ``"sn_wavefront"``
        and ``shape`` is the whole-interior flat shape
        ``(layout.total_size,)``. Inner product is Euclidean
        (``inner_product_weights=None``) — interior sweep-transient state.
    layout : FaceLayout
        The flat-buffer descriptor (one face-normal field per active axis;
        per-axis shapes; offsets). ``compare=False`` leaf-data so it does
        NOT pollute the ``(name, shape)`` identity — two interior face
        spaces of the same total interior size compare equal regardless of
        their per-axis decomposition.
    """

    layout: Optional[FaceLayout] = field(
        default=None, repr=False, compare=False,
    )

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    #
    # Explicit delegation (matching TensorProductSpace / DualSpace): the
    # @dataclass(frozen=True) decorator would otherwise auto-generate an
    # __eq__ that also compares ``inner_product_weights`` (an ndarray when
    # set → "truth value ambiguous"). Delegation restores the (name, shape)
    # identity convention; ``layout`` is already compare=False.

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    @property
    def face_names(self) -> tuple[str, ...]:
        r"""Ordered interior face-field names (one per active axis)."""
        if self.layout is None:
            return ()
        return tuple(self.layout.faces)

    # ── Construction ─────────────────────────────────────────────────

    @staticmethod
    def interior_layout(
        N: int, ng: int, dims: tuple[int, ...],  # noqa: N803 — N is the quad size
    ) -> FaceLayout:
        r"""Axis-parametric interior :class:`FaceLayout` for ``dims`` cells.

        One face-normal field per axis ``a``: shape ``(N, ng, *d)`` where
        ``d[b] = dims[b] + 1`` for ``b == a`` (the faces along ``a`` number
        one more than the cells) and ``d[b] = dims[b]`` otherwise. So a 2-D
        ``dims=(nx, ny)`` mesh gives ``x``:``(N, ng, nx+1, ny)`` +
        ``y``:``(N, ng, nx, ny+1)``; a 1-D ``dims=(nx,)`` gives
        ``x``:``(N, ng, nx+1)``; a 3-D ``dims=(nx, ny, nz)`` adds
        ``z``:``(N, ng, nx, ny, nz+1)``. The axis count is a parameter, so
        3-D is not a new code path (plan §3a′).
        """
        named: list[tuple[str, tuple[int, ...]]] = []
        for a in range(len(dims)):
            shape = [N, ng]
            for b in range(len(dims)):
                shape.append(dims[b] + 1 if b == a else dims[b])
            named.append((_AXIS_NAMES[a], tuple(shape)))
        return FaceLayout.from_named_shapes(named)

    @classmethod
    def from_mesh(cls, mesh: "SNMesh") -> "InteriorFaceSpace":
        r"""Build the interior face space sized to ``mesh``.

        Dispatches 1-D (``mesh.ny == 1`` → ``dims=(nx,)``) vs 2-D Cartesian
        (``dims=(nx, ny)``) on ``mesh.ny`` — the same signal the sweep /
        matvec use (``ny > 1`` ⇒ 2-D Cartesian). 3-D (``Mesh3D``) is not
        reachable here yet; the layout *builder* (:meth:`interior_layout`)
        already accepts a 3-element ``dims`` (the type is forward-built).
        """
        N = mesh.quad.N
        ng = mesh.ng
        dims = (mesh.nx,) if mesh.ny == 1 else (mesh.nx, mesh.ny)
        layout = cls.interior_layout(N, ng, dims)
        return cls(
            name="sn_wavefront",
            shape=(int(layout.total_size),),
            layout=layout,
        )

    def __repr__(self) -> str:
        return f"InteriorFaceSpace({self.name!r}, shape={self.shape})"

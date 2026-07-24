r"""Augmented geometry for S\ :sub:`N` discrete ordinates transport.

:class:`SNMesh` is axis-primary (C5.1, #225): its canonical spatial
representation is a tuple of :class:`~orpheus.transport.mesh.axis.Axis1D`, and it
precomputes the coordinate-specific streaming stencil used by the
transport sweep. Two construction surfaces funnel into one body — the
axis-native :meth:`SNMesh.from_axes`, and the legacy
:class:`~geometry.mesh.Mesh1D` / :class:`~geometry.mesh.Mesh2D`
constructor (converted to axes once at the boundary).

Three coordinate systems are supported: Cartesian (1D/2D), spherical
(1D), and cylindrical (1D).  Curvilinear geometries precompute angular
redistribution coefficients (:math:`\alpha`), the geometry factor
:math:`\Delta A/w`, and Morel--Montry angular closure weights.
"""

from __future__ import annotations

import warnings
from functools import cached_property
from typing import ClassVar, Iterator, TYPE_CHECKING

import numpy as np

from orpheus.geometry import CoordSystem, Mesh1D, Mesh2D
from orpheus.geometry.boundary import (
    BoundaryTraceLaw,
    ReflectiveBoundary,
    VacuumInflow,
)
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.geometry.reduced_operator import (
    ReducedStreamingOperator,
    StreamingTerms,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.transport.method import resolve_boundary_conditions
from orpheus.transport.mesh.axis import (
    Axis1D,
    FaceLabel,
    axes_from_legacy_mesh,
    face_labels as _axis_face_labels,
    face_outflow_ordinates as _axis_face_outflow_ordinates,
    face_shape as _axis_face_shape,
    legacy_mesh_from_axes,
    n_unknowns_flat as _axis_n_unknowns_flat,
)
from orpheus.transport.mesh.material_mesh import (
    InconsistentMaterialsError,
    MaterialMesh,
)
from ..boundary.realizer import SNBoundaryRealizer
from .method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial.scheme import DiscretizationSchemeBase, CellVisit
from orpheus.transport.spatial.diamond import DiamondDifference
from ..sweep.pole_angular_closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    PoleAngularClosureBase,
    default_angular_closure_class,
    morel_montry_tau_raw_per_level,
)

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.numerics.spaces.radial_characteristic_space import (
        RadialCharacteristicBoundarySpace,
        RadialCharacteristicInteriorSpace,
    )
    # NOTE (B.5.A): the mesh provides shape data only and does NOT import
    # transport-field types — zero-allocation lives on the field types
    # (``Field.zeros`` / ``<Leaf>.zeros_on`` / ``TimedFullField.zeros(...)``).
    # The ``AngularBoundaryFlux`` / ``AngularFlux`` mentions below are docstring
    # cross-references (Sphinx resolves them by full path, no import needed).


# ``InconsistentMaterialsError`` moved to
# :mod:`orpheus.transport.mesh.material_mesh` (it is raised by
# ``MaterialMesh.ng``, the method-agnostic group-consistency check) and is
# re-exported here for the SN-side consumers / tests that import it from
# ``orpheus.sn.mesh.augmented_mesh``.


# ═══════════════════════════════════════════════════════════════════════
# SNMesh
# ═══════════════════════════════════════════════════════════════════════

class SNMesh(MaterialMesh):
    """Augmented geometry for the discrete ordinates method.

    Axis-primary (C5.1, #225): the canonical spatial representation is
    :attr:`axes` — a tuple of :class:`~orpheus.transport.mesh.axis.Axis1D` — from
    which all shape metadata derives. Constructed either axis-natively
    via :meth:`from_axes` or from a legacy
    :class:`~geometry.mesh.Mesh1D` / :class:`~geometry.mesh.Mesh2D`
    (converted to axes once at the inbound boundary; the legacy object
    is retained as :attr:`mesh` for the consumers still reading through
    it). Precomputes the streaming stencil (diamond-difference
    coefficients that depend only on geometry and angular quadrature,
    not on cross sections).

    For Cartesian geometry the stencil stores one per-axis array, read via
    :meth:`streaming`:

    * ``streaming(a)[n, i]`` = :math:`2|\\mu_{a,n}| / \\Delta a_i`
      for every axis ``a < ndim`` (built over ``range(ndim)`` from
      ``quad.axis_cosines(a)`` — no hand-listed x/y pair).

    For future curvilinear geometries, additional curvature terms
    (:math:`\\alpha_n / r_i`) will be stored in ``self.curvature``.

    Parameters
    ----------
    mesh : Mesh1D or Mesh2D
        Base geometry.
    quadrature : Quadrature
        Angular quadrature (Gauss–Legendre, Lebedev, etc.).
    materials : dict mapping material id to Mixture
        Macroscopic cross sections keyed by the integer ids appearing
        in ``mesh.mat_ids`` / ``mesh.mat_map``.  Required (Issue #197
        PR-TYPED-0).  The authoritative source of truth for both
        cross sections and the group count :attr:`ng`; every operator
        that consumes ``sn_mesh`` (L, C, S, F) reads materials from
        here, not from a parallel argument.  All materials must agree
        on ``ng`` — heterogeneous group structures are a
        homogenization-step concern that must precede SNMesh
        construction.

    Attributes
    ----------
    materials : dict mapping material id to Mixture
        The materials dict passed at construction (single source of
        truth).
    ng : int
        Number of energy groups, derived from materials and validated
        for consistency.
    BOUNDARY_OPERATOR_REGISTRY : dict[str, type[BoundaryTraceLaw]]
        Supported boundary-condition kinds (Wave 8 / C8.3) — the SN
        law-admission table read by the shared
        :func:`~orpheus.transport.method.resolve_boundary_conditions`
        body (#290 P7b). Values are :class:`BoundaryTraceLaw`
        subclasses (``VacuumInflow``, ``ReflectiveBoundary``) realized
        per face by :meth:`realize_boundary_law` via
        :class:`SNBoundaryRealizer` for every supported mesh
        (1-D Cartesian, 1-D spherical, 1-D cylindrical, 2-D
        Cartesian) and wrapped in :class:`_BoundBoundaryOperator`
        for compatibility with the SN-side call surface.
    bc : dict[str, _BoundBoundaryOperator]
        Resolved BC operator per boundary face, keyed by the face
        name — the SAME keys as :attr:`boundary_face_layout` /
        ``angular_trace.layout.faces``, both derived from :attr:`face_labels`
        through the single-sourced
        :attr:`~orpheus.transport.mesh.axis.FaceLabel.face_name` crosswalk (C4,
        #220). Each value is a :class:`_BoundBoundaryOperator` shim
        wrapping the realized 1-arg :class:`LinearOperator` and
        carrying a ``kind`` tag so ``sn_mesh.bc["xmin"] == "vacuum"``
        style comparisons work. The face inventory IS the BC
        inventory: slab ``{"xmin", "xmax"}``; **a solid sphere /
        cylinder has only ONE entry** (``"xmax"``, the outer radius —
        the pole r=0 is the angular closure's regularity condition,
        not a BC face, so it has NO entry rather than a ``None``);
        2-D Cartesian all four faces.
    """

    BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]] = {
        "vacuum": VacuumInflow,
        "reflective": ReflectiveBoundary,
    }
    # Values are the LAW CLASSES themselves (not factory functions), looked up
    # by the shared TransportMethod resolve body (#290 P7b —
    # ``resolve_boundary_conditions`` owns the face loop and the tag → law
    # parse; ``realize_boundary_law`` below dispatches
    # :class:`SNBoundaryRealizer`), applied uniformly for 1-D Cartesian, 1-D
    # spherical, 1-D cylindrical, and 2-D Cartesian meshes.
    #
    # The 5 other kinds the realizer handles today (``white``, ``periodic``,
    # ``albedo``, ``prescribed_inflow``, ``mixed``) are NOT registered here —
    # adding them requires SN-sweep-side wiring (sweep cycles for periodic,
    # etc.).  Future expansion is mechanical: add the law class as a value
    # here, ensure the realizer dispatch handles it (it does), and add an
    # SN-side test that the sweep behaves correctly.

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        quadrature: Quadrature,
        materials: "dict[int, Mixture]",
        scheme: DiscretizationSchemeBase | None = None,
        pole_angular_closure: "type[PoleAngularClosureBase] | None" = None,
    ) -> None:
        # The legacy inbound surface: convert the Mesh1D / Mesh2D declaration
        # to the canonical axis tuple ONCE at the boundary, extract the one
        # payload the axes cannot carry (the material assignment — named
        # ``mat_ids`` on Mesh1D, ``mat_map`` on Mesh2D), and run the same
        # construction body as :meth:`from_axes`. Everything downstream derives
        # from ``self.axes``; ``self.mesh`` survives as inbound provenance for
        # the consumers still reading through it (1-D reduced streaming
        # construction, trace build, realizer metadata, MMS helpers).
        self._init_core(
            axes=axes_from_legacy_mesh(mesh),
            mesh=mesh,
            mat_map=mesh.mat_ids if isinstance(mesh, Mesh1D) else mesh.mat_map,
            quadrature=quadrature,
            materials=materials,
            scheme=scheme,
            pole_angular_closure=pole_angular_closure,
        )

    def _init_core(
        self,
        *,
        axes: tuple[Axis1D, ...],
        mesh: Mesh1D | Mesh2D | None,
        mat_map: np.ndarray | None,
        quadrature: Quadrature,
        materials: "dict[int, Mixture]",
        scheme: DiscretizationSchemeBase | None,
        pole_angular_closure: "type[PoleAngularClosureBase] | None",
    ) -> None:
        # The ONE construction body both surfaces funnel into (C5.1).
        #
        # ── Method-agnostic DATA block → MaterialMesh base ──
        # :meth:`MaterialMesh._init_data` sets ``self.mesh`` / ``self.materials``
        # / ``self.axes`` / ``self.axis_widths`` / ``self.mat_map`` /
        # ``self._volumes`` / ``self._areas`` / ``self.nx`` / ``self.coord`` and
        # runs the materials-consistency validation.  ``materials`` is REQUIRED:
        # SNMesh IS the SN phase space (mesh × quadrature × material group
        # structure); without materials ``.ng`` is undefined (Pattern 4 —
        # illegal states unrepresentable).
        MaterialMesh._init_data(
            self,
            axes=axes,
            mesh=mesh,
            mat_map=mat_map,
            materials=materials,
        )

        # ── SN method layer (BEHAVIOR atop the MaterialMesh data) ──
        self.quad = quadrature
        # Cell-update strategy. Defaults to :class:`DiamondDifference`, which
        # reproduces the inlined sweep math bit-identically (every regression
        # snapshot at ``tests/sn/regression/snapshots/`` was generated with DD
        # and matches bit-for-bit through ``self.scheme.update(...)``).  Pass
        # ``scheme=LinearDiscontinuous()`` etc. to select another.
        self.scheme: DiscretizationSchemeBase = (
            scheme if scheme is not None else DiamondDifference()
        )
        # Angular-redistribution closure.  The default is
        # :class:`MorelMontryAngularSweep` for curvilinear (the canonical
        # Hébert §3.9.4 per-cell M-M weighted-DD angular recurrence with the
        # Carlson coupled-pole seed) and :class:`IdentityAngularClosure` for
        # Cartesian (flat geometry has no §3.9.4 term).  Derivation + the
        # ERR-026 closure: curvilinear_numerics.rst
        # §sn-phase-d-carlson-coupled-pole-sweep.
        #
        # Instantiation is DEFERRED until after the coord dispatch below
        # populates ``self.reduced`` / ``self._volumes`` / ``self.axis_widths``
        # (the data the strategies bind to) — see the ``self.pole_angular_closure
        # = …`` line after the BC resolution.  The override is a CLASS, not an
        # instance: a closure binds to its mesh at construction (``cls(self)``),
        # and this mesh does not exist yet when the caller assembles the
        # constructor arguments (Pattern 4 — an unbound / foreign-bound closure
        # is now unspellable).
        self._user_supplied_closure = pole_angular_closure

        # (``self.axes`` / ``self.axis_widths`` / ``self.mat_map`` /
        # ``self._volumes`` / ``self._areas`` / ``self.nx`` /
        # ``self.coord`` and the materials-consistency validation are all
        # set by the ``MaterialMesh._init_data`` call above — the
        # method-agnostic data block.)

        # Dispatch stencil setup by coordinate system.
        #
        # Curvilinear connection-coefficient math (sphere / cylinder) lives
        # in :mod:`orpheus.geometry.reduced_operator` (Wave B Issue 6) so
        # MoC and CP can consume the same primitive — Cardinal Rule 2
        # forbids duplicating it on each solver-side mesh class.  The
        # Cartesian per-axis streaming stencils are SN-specific (DD
        # denominator precomputation) and stay local to ``_setup_cartesian``.
        #
        # ``self.reduced`` is the canonical accessor every downstream
        # consumer should bind to: ``sn_mesh.reduced.streaming_terms(
        # cell_idx, dir_idx, mu_level_idx)`` returns the per-(cell,
        # direction) packet a sweep cell update needs (the deprecated
        # ``@property`` accessors below still preserve the legacy names).
        # ``self.coord`` was derived from the axes by
        # ``MaterialMesh._init_data`` (the whole-mesh coordinate system;
        # multi-axis tuples are all-Cartesian by construction). The 1-D
        # arms hand the legacy ``Mesh1D`` to the reduced streaming
        # constructors (the genuine remaining Mesh1D consumers — shared
        # with MoC/CP via :mod:`orpheus.geometry.reduced_operator`).
        match self.coord:
            case CoordSystem.CARTESIAN:
                self._setup_cartesian()
                # Slab also gets a ``ReducedStreamingOperator`` for
                # completeness so ``sn_mesh.reduced`` is always populated;
                # the slab variant carries empty curvature arrays and
                # ``requires_upstream_angular_state = False``.
                if self.ndim == 1:
                    assert isinstance(mesh, Mesh1D)
                    self.reduced: ReducedStreamingOperator = slab_streaming(
                        mesh, quadrature,
                    )
                else:
                    self.reduced = None  # type: ignore[assignment]
            case CoordSystem.CYLINDRICAL:
                assert isinstance(mesh, Mesh1D)
                self.reduced = cylindrical_streaming(mesh, quadrature)
                self.curvature: str | None = "cylindrical"
                # Cartesian-style per-axis streaming arrays not used here
                # (curvilinear streaming lives in reduced.streaming_terms).
                self._streaming_axes = None
            case CoordSystem.SPHERICAL:
                assert isinstance(mesh, Mesh1D)
                self.reduced = spherical_streaming(mesh, quadrature)
                self.curvature = "spherical"
                self._streaming_axes = None

        # ── Boundary trace + realized laws ──
        # Build ONE unified trace space per SNMesh, keyed on the mesh's
        # TRUE boundary faces (``boundary_face_layout``): slab
        # ``xmin``/``xmax``, curvilinear ``xmax`` only (the pole at r=0
        # is the angular closure's regularity condition, not a BC
        # face), multi-D Cartesian all ``2·ndim`` faces. Inflow /
        # outflow are selectors over the signed Ω·n it carries.
        # UNCONDITIONAL — every constructible SNMesh builds its trace
        # (geometry-blind: quadrature + face names); built HERE, in the
        # construction body, as phase-space substrate — not inside BC
        # resolution.
        from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
        self._trace = AngularTraceSpace.from_quadrature_and_layout(
            self.quad, self.boundary_face_layout,
        )

        # Resolve the per-axis BC declarations through the ONE shared
        # TransportMethod body (#290 P7b): the face loop over
        # ``face_labels``, the ``BC("reflective")`` infinite-lattice /
        # eigenvalue default, and the tag → law parse are method-
        # generic; :meth:`realize_boundary_law` below is the SN arm.
        # The face inventory IS the BC inventory by construction (C4,
        # #220): a face that exists has exactly one entry; the curvilinear
        # pole has none (Pattern 4 — a pole-BC is unrepresentable).
        # Consumers key into :attr:`bc` by the same face names the trace
        # layout carries.
        self.bc: dict[str, _BoundBoundaryOperator] = (
            resolve_boundary_conditions(self)
        )

        # (Materials-consistency validation — every ``mat_map`` id has a
        # ``materials`` entry, and all materials agree on ``ng`` — plus
        # the eager ``ng`` trigger run inside ``MaterialMesh._init_data``
        # above, so a bad materials dict raises at construction time.)

        # ── Pole-angular closure binding (PR-TYPED-6.5 Phase 2.9) ──
        # All upstream state needed by the closure constructors is now
        # available (``self.reduced``, ``self._volumes``, ``self.axis_widths``,
        # ``self.quad``, ``self.ng``).  The user-supplied closure CLASS —
        # or the default-by-coord-system — binds to ``self`` through the
        # family's one-positional-mesh construction contract
        # (``cls(sn_mesh)``); every mesh therefore carries a BOUND closure.
        closure_cls = (
            self._user_supplied_closure
            if self._user_supplied_closure is not None
            else default_angular_closure_class(self.coord)
        )
        self.pole_angular_closure: PoleAngularClosureBase = closure_cls(self)
        # Drop the temporary attribute now that the closure is bound.
        del self._user_supplied_closure

    # ── Boundary condition resolution ─────────────────────────────────
    #
    # The face loop, the reflective default, and the tag → law parse
    # live in the ONE shared TransportMethod body,
    # :func:`~orpheus.transport.method.resolve_boundary_conditions`
    # (#290 P7b — it replaced the twin ``SNMesh._resolve_bcs`` /
    # ``DiffusionMesh._resolve_bcs`` loops). Only the genuinely
    # SN-specific arm remains here:

    def realize_boundary_law(
        self,
        law: BoundaryTraceLaw,
        face: str,
    ) -> "_BoundBoundaryOperator":
        r"""Realize one typed boundary law on ``face`` — the SN arm of the
        :class:`~orpheus.transport.method.TransportMethod` hook.

        Called per face by the shared
        :func:`~orpheus.transport.method.resolve_boundary_conditions`
        body. Build an :class:`SNMethodSpace` carrying the precomputed
        unified :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`,
        hand the law to :class:`SNBoundaryRealizer.realize`, wrap the
        1-arg result in :class:`_BoundBoundaryOperator` so the SN-side
        call surface sees a uniform 1-arg ``apply(psi)`` contract with
        the ``kind`` string tag.

        The ``kind`` tag reads the LAW's own registry key
        (:attr:`~orpheus.numerics.registry.RegistryMixin.key` —
        ``"vacuum"``, ``"reflective"``): the kind string's single
        source is the law class itself, and every
        ``BOUNDARY_OPERATOR_REGISTRY`` entry maps a tag to the law
        registered under that same key, so the tag equals the declared
        ``BC.kind`` by construction.

        Issue #188 / C188.3: every supported mesh (1-D Cartesian,
        1-D spherical, 1-D cylindrical, 2-D Cartesian) routes
        through the realizer here. The pre-C188.3 curvilinear
        bypass — which wrapped the raw 2-arg
        :class:`BoundaryTraceLaw` with a bound quadrature — is
        gone, made redundant by the unified trace's curvilinear
        support. ``face`` must name a face present in the trace;
        curvilinear's inner pole has no label and is handled by the
        angular closure, not here.
        """
        method_space = SNMethodSpace.for_face(
            mesh=self.mesh,
            quadrature=self.quad,
            face=face,
            trace=self._trace,
        )
        realized = SNBoundaryRealizer().realize(law, method_space)
        return _BoundBoundaryOperator(realized, kind=law.key)

    # ── Properties ────────────────────────────────────────────────────
    #
    # ``_validate_materials`` and the data properties ``ng`` / ``volumes``
    # / ``volume_measure`` / ``areas`` / ``ndim`` / ``spatial_shape`` —
    # plus the ``material_xs_field()`` builder — are inherited from
    # :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` (the
    # method-agnostic data carrier).  SNMesh adds only the SN-method
    # behavior (quadrature / streaming stencil / boundary trace / closures)
    # on top.

    @property
    def is_1d(self) -> bool:
        """True if this is a genuine 1-D mesh (``ndim == 1``).

        Reads the genuine spatial dimensionality, NOT the phantom
        ``ny == 1`` shim: a :class:`Mesh2D` with a single y-cell is 2-D
        (``ndim == 2``) and is NOT 1-D. The old ``ny == 1`` test
        misclassified that degenerate case and was the root of #214; the
        genuine-dimensionality test is the phantom-axis-elimination
        invariant (R-1 Phase A). This is the single source of truth for
        the 1-D-vs-multi-D dispatch in the streaming operators
        (``not sn_mesh.is_1d`` selects the multi-D Cartesian path).
        """
        return self.ndim == 1

    @property
    def is_cartesian(self) -> bool:
        """True if the mesh carries no curvature (Cartesian slab / 2-D / 3-D).

        The genuine coordinate-system criterion — ``curvature is None`` for a
        Cartesian slab or a multi-D Cartesian mesh; a string
        (``'spherical'`` / ``'cylindrical'``) for a curvilinear 1-D mesh.
        This is ORTHOGONAL to :attr:`is_1d`: a slab is Cartesian AND 1-D; a
        2-D Cartesian mesh is Cartesian AND not 1-D; a cylinder is 1-D AND
        not Cartesian.  Sweep-strategy selection
        (``orpheus.sn.loss_representation.default_for``) keys on BOTH axes —
        the anti-hyperplane DAG family requires ``is_cartesian``, the chain
        scan requires ``is_1d`` — so neither alone is a sufficient
        discriminator.
        """
        return self.curvature is None

    def streaming(self, axis: int) -> np.ndarray:
        r"""Per-axis RAW down-face streaming ``g = |μ_axis|·A_down/V = |μ_axis|/Δ_axis``, ``(N, n_axis)``.

        The dimension-generic accessor the anti-hyperplane DAG walk reads as
        ``str_axes[axis]`` — the **scheme-agnostic** geometric streaming.  Each
        spatial scheme applies its OWN closure factor: DD contributes
        :math:`\sum_a 2g_a` to the cell-balance denominator (the
        :math:`2 = 1/w_{\rm DD}` is DD's diamond closure, applied in its kernel,
        NOT here — #240); Linear-Discontinuous reads the same raw ``g``.
        Indexes the per-axis stencil tuple ``_setup_cartesian`` builds over
        ``range(ndim)`` (since C3.6 there is no hand-listed
        ``(streaming_x, streaming_y)`` pair to drift out of axis order — the
        tuple IS positional-by-axis from birth).

        Cartesian-only (the anti-hyperplane lattice is a Cartesian object);
        curvilinear meshes carry their streaming in
        ``reduced.streaming_terms`` (the chain-scan substrate) and are swept by
        the ``CumprodScan`` strategy, not the DAG walk.
        ``axis`` must satisfy ``0 <= axis < ndim``.
        """
        # The None-ness IS the Cartesian-only gate: ``_setup_cartesian``
        # builds the stencil tuple, the curvilinear arms assign ``None`` —
        # so checking the attribute directly both guards and narrows.
        streaming_axes = self._streaming_axes
        if streaming_axes is None:
            raise AttributeError(
                "SNMesh.streaming(axis) is Cartesian-only; curvilinear meshes "
                "carry streaming in reduced.streaming_terms (the chain-scan "
                "substrate), not the anti-hyperplane DAG."
            )
        if not 0 <= axis < self.ndim:
            raise IndexError(
                f"streaming axis {axis} out of range for ndim={self.ndim}"
            )
        return streaming_axes[axis]

    # ── Dim-agnostic geometry primitives (R-1 Phase A C1) ─────────────
    #
    # ``ndim`` / ``spatial_shape`` are inherited from
    # :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` (the
    # method-agnostic data carrier).

    @property
    def face_labels(self) -> tuple[FaceLabel, ...]:
        r"""Canonical boundary-face inventory.

        Each :class:`~orpheus.transport.mesh.axis.FaceLabel` carries an
        ``axis_index`` and an ``endpoint`` label, derived from the
        per-axis endpoints. Cartesian 1-D returns 2 labels; spherical
        / cylindrical 1-D returns 1 label (the pole is NOT a face —
        see :class:`~orpheus.transport.mesh.axis.Axis1D` docstring); 2-D Cartesian
        returns 4 labels; synthetic 3-D Cartesian would return 6.

        The iteration order is the canonical concatenation order for
        :meth:`AngularFlux.to_flat` (C3) and the canonical iteration
        order for :attr:`AngularBoundaryFlux.face_buffers` (C4).
        """
        return _axis_face_labels(self.axes)

    def face_shape(self, label: FaceLabel) -> tuple[int, ...]:
        r"""Spatial shape of the boundary face identified by ``label``.

        The face lies in the codimension-1 hyperplane spanned by the
        axes other than ``label.axis_index``; its shape is the
        per-axis cell count of those axes in axis-index order.
        """
        return _axis_face_shape(self.axes, label)

    def face_outflow_ordinates(self, label: FaceLabel) -> np.ndarray:
        r"""Ordinate indices whose direction-cosine is OUTWARD at face ``label``.

        At the ``max`` / ``outer`` endpoint of an axis, outflow is
        :math:`\mu_{axis} > +10^{-15}`; at ``min``, outflow is
        :math:`\mu_{axis} < -10^{-15}`. Ordinates exactly tangent to
        the face contribute neither inflow nor outflow.

        This method is the canonical producer for the per-face
        outflow mask used by the pack convention (C3),
        :class:`AngularBoundaryFlux.face_buffers` (C4), and the sweep DAG
        face-trace state (C5).
        """
        return _axis_face_outflow_ordinates(self.axes, label, self.quad)

    @property
    def n_unknowns_flat(self) -> int:
        r"""Total flat-vector size for typed :class:`AngularFlux`.

        The pack convention (C3) is the direct-sum decomposition
        :math:`V = V_\text{cells} \oplus \bigoplus_\ell V_{\text{face}, \ell}`;
        ``n_unknowns_flat`` is the dimension of that vector space.
        Cells contribute :math:`N \cdot n_g \cdot \prod_i n_i`; each
        face :math:`\ell` contributes
        :math:`|\text{outflow}_\ell| \cdot n_g \cdot \prod_{i \ne \text{axis}(\ell)} n_i`.
        """
        return _axis_n_unknowns_flat(self.axes, self.quad, self.ng)

    @classmethod
    def from_axes(
        cls,
        axes: tuple[Axis1D, ...],
        quadrature: "Quadrature",
        materials: "dict[int, Mixture]",
        *,
        mat_map: np.ndarray | None = None,
        scheme: DiscretizationSchemeBase | None = None,
        pole_angular_closure: "type[PoleAngularClosureBase] | None" = None,
    ) -> "SNMesh":
        r"""Build an :class:`SNMesh` from an axis tuple — the axis-native surface.

        C5.1 (axis-primary inversion, #225): the caller's axes ARE the
        mesh's axes — stored verbatim and never round-tripped through a
        legacy mesh and re-derived (the pre-C5.1 round-trip silently
        reset custom endpoint labels to ``min``/``max``/``outer``). A
        legacy :class:`Mesh1D` / :class:`Mesh2D` ADAPTER is still
        synthesized at d≤2 for the consumers that read through
        ``self.mesh`` (1-D reduced streaming construction, trace build,
        realizer metadata) — each dissolves across C5.2–C5.5.

        Endpoint labels must be canonical (``min``/``max``/``outer``):
        the :attr:`bc` dict is keyed by
        :attr:`~orpheus.transport.mesh.axis.FaceLabel.face_name`, which fails loud
        on a custom label (C4 doctrine — overridable labels cannot
        silently desync the face-name crosswalk). Custom labels are for
        standalone axis use, not SNMesh construction.

        Parameters
        ----------
        axes : tuple of :class:`~orpheus.transport.mesh.axis.Axis1D`
            Per-axis 1-D mesh descriptors. Length 1 → 1-D mesh;
            length 2 → 2-D Cartesian mesh; length ≥3 → d-D Cartesian
            (C5.5, #225 — all-Cartesian required, mesh-adapter-free
            from birth, swept by the d-generic ``FullFieldWavefront``
            spine).
        quadrature : :class:`Quadrature`
            Angular quadrature.
        materials : dict[int, Mixture]
            Materials dict keyed by material id; same contract as the
            legacy constructor.
        mat_map : ndarray or None
            Material-id assignment. Shape ``spatial_shape``. Defaults
            to all-zeros (single material with id 0).
        scheme : DiscretizationSchemeBase or None
            Cell-update strategy. Defaults to :class:`DiamondDifference`.
        pole_angular_closure : type[PoleAngularClosureBase] or None
            Override the default pole-angular closure CLASS
            (curvilinear → :class:`MorelMontryAngularSweep`,
            Cartesian → :class:`IdentityAngularClosure`).  A class, not
            an instance: closures bind to their mesh at construction
            (``cls(sn_mesh)``), and the mesh does not exist yet.
        """
        axes = tuple(axes)
        # C5.5 (#225): d≥3 is mesh-adapter-free from birth — every
        # consumer that read through ``self.mesh`` was dissolved across
        # C5.2–C5.4 (volume measure, trace, windowing gates) or is
        # d≤2-only (the 1-D reduced streaming constructors, the MMS
        # helpers). d≤2 still synthesizes the legacy adapter for those
        # remaining consumers.
        mesh = (
            legacy_mesh_from_axes(axes, mat_map=mat_map)
            if len(axes) <= 2 else None
        )
        obj = cls.__new__(cls)
        obj._init_core(
            axes=axes,
            mesh=mesh,
            mat_map=mat_map,
            quadrature=quadrature,
            materials=materials,
            scheme=scheme,
            pole_angular_closure=pole_angular_closure,
        )
        return obj

    @classmethod
    def from_material_mesh(
        cls,
        material_mesh: MaterialMesh,
        quadrature: "Quadrature",
        *,
        scheme: DiscretizationSchemeBase | None = None,
        pole_angular_closure: "type[PoleAngularClosureBase] | None" = None,
    ) -> "SNMesh":
        r"""Promote a :class:`MaterialMesh` to a solvable SN phase space.

        The data/behavior join: a :class:`MaterialMesh` carries the
        method-agnostic data (axes + materials + mat_map + volumes); this
        classmethod adds the SN method layer (angular quadrature + sweep
        stencil + boundary trace + closures) to make it solvable.

        It is the natural consumer of cross-section homogenization: a
        homogenized :class:`MaterialMesh` (from
        :meth:`~orpheus.sn.solution.Solution.homogenize`) is promoted
        here to re-solve the coarsened problem on the same outer geometry
        (the "re-solve the homogenized problem" path).  The
        material-mesh's axes / mesh / mat_map / materials are passed
        through verbatim — ``_init_core`` re-derives the data block
        bit-identically from them.

        Parameters
        ----------
        material_mesh : MaterialMesh
            The mesh+materials carrier to promote.
        quadrature : Quadrature
            Angular quadrature for the SN method.
        scheme : DiscretizationSchemeBase or None
            Cell-update strategy.  Defaults to :class:`DiamondDifference`.
        pole_angular_closure : type[PoleAngularClosureBase] or None
            Override the default pole-angular closure CLASS
            (curvilinear → :class:`MorelMontryAngularSweep`,
            Cartesian → :class:`IdentityAngularClosure`).  A class, not
            an instance: closures bind to their mesh at construction
            (``cls(sn_mesh)``), and the mesh does not exist yet.
        """
        obj = cls.__new__(cls)
        obj._init_core(
            axes=material_mesh.axes,
            mesh=material_mesh.mesh,
            mat_map=material_mesh.mat_map,
            quadrature=quadrature,
            materials=material_mesh.materials,
            scheme=scheme,
            pole_angular_closure=pole_angular_closure,
        )
        return obj

    @property
    def angular_trace(self) -> "AngularTraceSpace":
        r"""The unified boundary :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`.

        One concrete trace space for the whole boundary :math:`\Gamma`,
        built (A.2/A.3) from this mesh's quadrature +
        :attr:`boundary_face_layout` (geometry-blind since C5.3, #225).
        It is the single source of truth for the signed projection
        :math:`\Omega\cdot\hat n_f` per face; the inflow / outflow
        *selectors*
        (:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`)
        replace the inline ``sign(Ω·n)`` masks that the streaming matvec
        and the boundary realizer previously each recomputed.

        ALWAYS non-``None`` (C5.3): the only mesh the pre-C5.3 gate
        excluded — a cylindrical :class:`~orpheus.geometry.mesh.Mesh2D`
        — cannot become an SNMesh at all, so every constructible SNMesh
        carries a trace.
        """
        return self._trace

    #: FP-noise guard for the R12a strict-interior test on the first-ordinate
    #: raw M-M weight. The production instances are BIT-EXACT members of
    #: {0, 1} (product rules: node ON the starting edge → τ_raw = 0;
    #: level-symmetric rules: duplicate-η midpoint edge collapses onto η₀ →
    #: τ_raw = 1) or safely interior (sphere-GL ≈ 0.39–0.42), so the guard
    #: never decides a real case — it exists so a future quadrature whose
    #: coincidence arithmetic differs by ULPs cannot flip presence, mirroring
    #: the closure's own ``abs(dμ) > 1e-15`` degenerate-cell guard.
    _SEED_TAU_EPS: ClassVar[float] = 1e-12

    @cached_property
    def radial_characteristic_levels(self) -> tuple[int, ...]:
        r"""μ-level indices that consume INDEPENDENT starting-direction state (R12a).

        The seed-presence predicate of #282 route (a): a level carries a
        ψ½ block iff its first-ordinate raw Morel–Montry weight
        :math:`\tau_{\rm raw,0}` lies strictly in :math:`(0, 1)` — i.e.
        the M-M half-angle recurrence genuinely consumes a seed value
        that is neither a rank-duplicate of the level's first node
        (:math:`\tau_{\rm raw} = 0`, cylinder product rules — the #229
        clamp fact) nor dead under the recurrence's :math:`(1-\tau_0)`
        thread weight (:math:`\tau_{\rm raw} = 1`, cylinder
        level-symmetric rules — duplicate-η midpoint edges). Sphere-GL
        is the carrying instance (one level, the whole quadrature,
        :math:`\tau_{\rm raw,0} \approx 0.39\text{–}0.42`); Cartesian
        never carries.

        R12a refines the R12 letter ("μ_start ∉ the level's μ-nodes"),
        whose claimed ⟺ with ``τ_raw ≠ 0`` fails on level-symmetric
        cylinder rules (μ_start ∉ nodes there, yet the seed is dead —
        measured 0.0-bit solve insensitivity). Single-sourced from
        :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_raw_per_level`
        — the SAME edge construction the production τ clamp consumes.
        Level indexing matches the closure's: the sphere's single M-M
        level is index ``0``; cylinder levels index
        ``quad.level_indices``.
        """
        if self.curvature is None:
            return ()
        assert self.reduced is not None  # curvilinear ⇒ reduced populated
        raw = morel_montry_tau_raw_per_level(self.quad, self.reduced.coord)
        eps = self._SEED_TAU_EPS
        return tuple(
            p for p, tau_level in enumerate(raw)
            if eps < float(tau_level[0]) < 1.0 - eps
        )

    def _radial_characteristic_for_levels_args(
        self,
    ) -> "tuple[tuple[int, ...], int, int, np.ndarray] | None":
        r"""The shared ``for_levels`` args for the ψ½ ray spaces, or ``None``.

        Single-sources the R12a levels gate + the ``(ng, nx, cell_volumes)``
        sourcing across the split :attr:`radial_characteristic_interior_space`
        / :attr:`radial_characteristic_boundary_space` (Phase B — the
        coupled-block campaign poses the ψ½ ray as System B, its own
        interior ⊕ boundary composite; the unified space retired at 4e), so
        the ψ½ spaces are built from ONE set of inputs. ``None`` ⟺
        :attr:`radial_characteristic_levels` is empty (absence is spelled
        ``None``, never a zero-DOF space). ``cell_volumes`` is the
        ``G_sd = V_cell`` state metric — the SAME radial cell-volume measure the
        bulk metric ``G_bulk = V_cell·w_n`` reads (:attr:`full_field_space`).
        """
        levels = self.radial_characteristic_levels
        if not levels:
            return None
        return (
            levels,
            self.ng,
            int(self.spatial_shape[0]),
            np.asarray(self.volumes, dtype=float).ravel(),
        )

    @cached_property
    def radial_characteristic_interior_space(
        self,
    ) -> "RadialCharacteristicInteriorSpace | None":
        r"""The ψ½ INTERIOR (cells) space — System B's interior block, or ``None``.

        The ``(ng, nx)`` cells legs under the SPD ``G_sd = V_cell`` metric, on
        which
        :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
        (A_BB) marches — the seed sibling of :attr:`angular_trace`, paired
        with :attr:`radial_characteristic_boundary_space` (Phase B; the
        historical unified space retired at 4e). ``None`` on non-carrying
        meshes (R12a). Cached (one per mesh), built from the shared
        :meth:`_radial_characteristic_for_levels_args`.
        """
        args = self._radial_characteristic_for_levels_args()
        if args is None:
            return None
        from orpheus.numerics.spaces.radial_characteristic_space import (
            RadialCharacteristicInteriorSpace,
        )

        levels, ng, nx, cell_volumes = args
        return RadialCharacteristicInteriorSpace.for_levels(
            levels, ng=ng, nx=nx, cell_volumes=cell_volumes,
        )

    @cached_property
    def radial_characteristic_boundary_space(
        self,
    ) -> "RadialCharacteristicBoundarySpace | None":
        r"""The ψ½ BOUNDARY (corner) space — System B's boundary block, or ``None``.

        The corner sibling of :attr:`radial_characteristic_interior_space`
        (Phase B): the
        ``(ng,)`` r = R corner legs under the ``G = V(r = R)`` corner gauge, on
        which
        :class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`
        (B_b) acts (inflow = given data; outflow = the defect row). ``None`` on
        the same non-carrying meshes (R12a). Cached (one per mesh), built from
        the shared :meth:`_radial_characteristic_for_levels_args`.
        """
        args = self._radial_characteristic_for_levels_args()
        if args is None:
            return None
        from orpheus.numerics.spaces.radial_characteristic_space import (
            RadialCharacteristicBoundarySpace,
        )

        levels, ng, nx, cell_volumes = args
        return RadialCharacteristicBoundarySpace.for_levels(
            levels, ng=ng, nx=nx, cell_volumes=cell_volumes,
        )

    @cached_property
    def radial_characteristic_field_space(self) -> "FullFieldSpace | None":
        r"""System B's member space — the ψ½ ``interior ⊕ boundary`` composite, or ``None``.

        The :class:`~orpheus.numerics.space.FunctionSpace` the re-typed
        coupling blocks declare (B.2b DP1): ``A_BA``'s codomain and ``B_b``'s
        domain/codomain — the carrier space of
        :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`.
        REUSES the family-blind
        :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
        (the same direct-sum member-wise metric dispatch System A's
        :attr:`full_field_space` uses — zero new space classes; this IS the
        post-eviction end-state, one composite-space class with instances
        differing in members), instantiated over the two split ψ½ spaces:

        * **interior** — :attr:`radial_characteristic_interior_space`
          (``G_sd = V_cell`` cells state metric);
        * **trace slot** — :attr:`radial_characteristic_boundary_space`
          (the ``G = V(r = R)`` corner gauge).

        Identity ``("radial_characteristic", (n_interior + n_corner,))`` —
        the name signals the instance. ``None`` on non-carrying meshes
        (R12a; System B does not exist there). Cached: one space per mesh,
        so every block shares one identity instance.
        """
        interior = self.radial_characteristic_interior_space
        boundary = self.radial_characteristic_boundary_space
        if interior is None or boundary is None:
            return None
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        return FullFieldSpace.from_blocks(
            interior, boundary, name="radial_characteristic",
        )

    @cached_property
    def full_field_space(self) -> "FullFieldSpace":
        r"""The composite carrier :math:`V_{\rm bulk} \oplus V_{\rm trace}` (Wave O / O.2b).

        The function space of the FULL streaming operator
        (:class:`~orpheus.sn.operators.streaming.StreamingOperator`) and every bulk
        :math:`\oplus` boundary composite — the domain/codomain under which
        ``L.H`` and ``(L + C - B).H`` are the **metric-correct G-adjoint**
        :math:`A^\dagger = G^{-1} A^{\mathsf T} G` (Issue #208). The
        block-diagonal Hilbert metric :math:`G` is

        * **bulk** :math:`G_{\rm bulk} = V_{\rm cell}\,w_n` — the full
          phase-space measure :math:`\mathrm{d}V\,\mathrm{d}\Omega`, stored
          ``(N, 1, nx, ny)`` so the per-ordinate angular weight ``w_n``
          (axis 0) and the per-cell volume ``V`` (axes 2,3) broadcast
          across the energy-group axis against the ``(N, ng, nx, ny)`` bulk
          tensor;
        * **trace** :math:`G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n` — the
          partial-current surface metric already carried by
          :attr:`angular_trace`.

        The two factors carry :math:`w_n`; they differ only in the
        spatial measure (cell volume vs. oriented face). On a carrying
        mesh the ψ½ ray's state metric (``G_sd = V_cell``) lives on
        **System B's own composite space**,
        :attr:`radial_characteristic_field_space` — never as a third
        block here (the B.2d eviction; the coupled DOF count is the
        honest two-system sum). Cached: the composite is immutable for a
        given mesh + quadrature.
        """
        from orpheus.numerics.space import FunctionSpace
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        N = self.quad.N
        w_n = np.asarray(self.quad.weights, dtype=float)  # (N,)
        V = np.asarray(self.volumes, dtype=float)         # (*spatial)
        # G_bulk = V_cell · w_n on (N, 1, *spatial) — group-independent
        # (broadcast over the ng axis): w_n along the ordinate axis, V along
        # the spatial axes (the full phase-space measure dV·dΩ). Rank-generic
        # over ndim so 1-D → (N, 1, nx), 2-D → (N, 1, nx, ny).
        g_bulk = (
            w_n.reshape((N, 1) + (1,) * V.ndim)
            * V.reshape((1, 1) + V.shape)
        )
        bulk_shape: tuple[int, ...] = (N, self.ng, *self.spatial_shape)
        # A multi-moment closure (LD) carries the trailing 2^d spatial-moment
        # axis on the bulk field, and its Hilbert metric MUST carry the
        # scheme's moment mass ``M_ii/V = ∏_a θ^{o_a}`` on that axis (#310 C2
        # ruling 3): ``G_bulk = V·w_n ⊗ diag(1, θ, …)``.  An average-only
        # ``V·w_n`` broadcast over the moment axis mis-weights the slope DOF
        # — ``.H`` becomes a WRONG adjoint on the slope rows AND reciprocity
        # goes Mode-12-blind to a slope-row transpose.  DD/Step
        # (``per_axis == 1``) take neither branch — shape and weights
        # byte-identical.
        if self.scheme.spatial_basis_per_axis > 1:
            moment_mass = self.scheme.moment_mass_diagonal(self.ndim)
            bulk_shape = bulk_shape + (moment_mass.size,)
            g_bulk = g_bulk[..., None] * moment_mass
        interior_space = FunctionSpace(
            name="sn_bulk",
            shape=bulk_shape,
            inner_product_weights=g_bulk,
        )
        return FullFieldSpace.from_blocks(
            interior_space,
            self.angular_trace,
        )

    @property
    def boundary_face_layout(self) -> "FaceLayout[str]":
        r"""Flat :class:`~orpheus.numerics.face_layout.FaceLayout` of boundary faces.

        Depth B step D-G primitive. Returns the per-geometry boundary
        face descriptor: which faces exist, their per-face shapes, and
        the flat-buffer offsets that pack them. The post-D-G pure-Field
        :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
        consumes this layout to lay out its flat backing buffer.

        Derived from :attr:`face_labels` (C4): one slot per label, named
        by the single-sourced :attr:`FaceLabel.face_name` crosswalk,
        shaped ``(N, ng, *face_shape(label))`` — axis-count generic, no
        per-geometry hand-list. The geometry mapping falls out:

        * 1-D slab — two faces ``xmin`` / ``xmax``, each ``(N, ng)``.
        * 1-D curvilinear sphere / cylinder — one face ``xmax``, shape
          ``(N, ng)`` (the single ``"outer"`` endpoint renders as
          ``xmax``; the geometric pole at r=0 is a regularity
          condition, not a BC face, so it has no label and no slot).
        * 2-D Cartesian — four faces: ``xmin`` / ``xmax`` shape
          ``(N, ng, ny)``; ``ymin`` / ``ymax`` shape ``(N, ng, nx)``.
        * A 3-axis mesh (C5) would yield six slots ``xmin`` … ``zmax``
          with codimension-1 shapes — no edit needed here.

        Spatial-moment tail (#251 — Leg B of #247)
        ------------------------------------------
        A multi-moment closure (LD's bilinear UBLD face) carries a
        trailing ``2^{d-1}``-transverse-moment axis per face slot, so a
        moment-resolved prescribed inflow can carry the along-face
        (transverse) Legendre slope and the sweep outflow can STORE
        those ``2^{d-1}`` moments instead of collapsing to the average
        (slot 0).  The width is the scheme's per-face count
        ``per_axis^{d-1}``, appended via the single-source
        :func:`~orpheus.numerics.moment_layout.face_moment_tail` (the
        same "append iff > 1" policy the cell-cochain
        :attr:`_LossRepresentation._n_face_moments` /
        :attr:`_spatial_moment_tail` key on).  DD/Step
        (``per_axis == 1`` → ``per_axis^{d-1} == 1`` → ``face_moment_tail
        == ()``) leaves every slot shape untouched, so the trace stays
        byte-identical — the negative control.  A 1-D slab face is a
        point (``face_shape == ()`` → ``per_axis^0 == 1`` → no tail even
        for LD-1D), so the transverse face-moment is a 2-D-and-higher
        concern by construction.

        Returns
        -------
        FaceLayout
            Per-geometry face descriptor. Total flat size = sum of
            ``prod(shape)`` over all faces. Slot order = the canonical
            :attr:`face_labels` order (axis ascending, endpoint in axis
            order), which reproduces the historical hand-listed order.

        Notes
        -----
        The layout contains ONLY boundary face slots. Interior
        wavefront cache cells (pre-D-G stored in AngularBoundaryFlux's 2-D
        ``xmin_xmax_buf`` / ``ymin_ymax_buf`` interior positions) are
        explicitly excluded — they live on
        :class:`~orpheus.sn.sweep_scratch.SweepScratch` post-D-G.
        """
        from orpheus.numerics.face_layout import FaceLayout
        from orpheus.numerics.moment_layout import face_moment_count, face_moment_tail

        N = self.quad.N
        # Per-face transverse moment count per_axis^{d-1} (#251) — the FACE
        # tail (the CELL tail is per_axis^d).  DD/Step → () → byte-identical.
        # Single-sourced with the cochain's ``_n_face_moments`` via
        # ``face_moment_count`` so the producer and the consumer cannot drift.
        n_face_moments = face_moment_count(self.scheme.spatial_basis_per_axis, self.ndim)
        moment_tail = face_moment_tail(n_face_moments)
        return FaceLayout.from_named_shapes([
            (label.face_name, (N, self.ng, *self.face_shape(label), *moment_tail))
            for label in self.face_labels
        ])

    # ── Sweep DAG traversal ───────────────────────────────────────────

    _DEGENERATE_ABS_ETA_THRESHOLD: ClassVar[float] = 1e-15

    def dag_walk(
        self,
        *,
        ordinate_idx: int | None = None,
        direction_sign: int | None = None,
        mu_level_idx: int | None = None,
    ) -> Iterator[CellVisit]:
        r"""Walk the per-ordinate cell DAG in topological order.

        Issue #196 Phase G Step 2.6 (Q3): the single canonical iteration
        primitive for 1-D sweeps.  Yields visits either for a single
        ordinate or for all ordinates of a sweep direction under one
        XOR signature.

        The SN sweep on a given ordinate is forward substitution on
        the block-triangular streaming + collision operator under the
        ordinate's DAG ordering.  This method yields the per-cell
        visit packets in that DAG order; the consumer folds over the
        packets, threading the spatial-upstream face flux through the
        accumulator and writing the per-cell angular state into a
        persistent array.

        Exactly one of ``ordinate_idx`` or ``direction_sign`` must be
        supplied (XOR):

        * ``ordinate_idx=n`` — yields visits for a single ordinate.
          For slab/sphere: ``n`` is the global ordinate index.  For
          cylindrical: ``n`` is the within-level azimuthal index
          :math:`m \in [0, M)` and ``mu_level_idx`` MUST also be
          supplied; the signed :math:`\eta` resolves through
          ``quad.level_indices[mu_level_idx][n]``.
        * ``direction_sign=±1`` — yields visits for the sweep
          direction (``+1`` outward, ``-1`` inward).  Cell ordering
          depends ONLY on the direction sign (and level for
          cylindrical), so any ordinate in the correct sign class
          yields the same cell sequence; this branch picks a
          non-degenerate representative.

        SN-specific by design.  MoC will not consume this method —
        its mathematical structure is fiber bundles + solution
        sheaves, a different DAG shape.  Premature abstraction
        avoided per Cardinal Rule 2.

        Parameters
        ----------
        ordinate_idx : int | None
            See above; mutually exclusive with ``direction_sign``.
        direction_sign : int | None
            See above; mutually exclusive with ``ordinate_idx``.
        mu_level_idx : int | None
            For cylindrical geometry: which :math:`\mu`-level the
            ordinate (subset) belongs to.  ``None`` for slab/sphere;
            required for cylindrical.

        Yields
        ------
        CellVisit
            One per cell, in topological order.  The packet's
            :attr:`face_area_downstream` is float: ``1.0`` for slab,
            ``0.0`` for cylindrical pure-azimuthal degenerate,
            physical face area for sphere / non-degenerate cylinder.

        Raises
        ------
        ValueError
            If neither or both of ``ordinate_idx`` / ``direction_sign``
            are supplied; if ``direction_sign not in (+1, -1)``; if
            called on a 2-D Cartesian mesh (no
            :class:`ReducedStreamingOperator`); if a cylindrical
            mesh is queried without ``mu_level_idx``; or if no
            non-degenerate representative ordinate exists for
            ``direction_sign``.

        Notes
        -----
        2-D Cartesian wavefront scheduling is intentionally not
        encapsulated here — its anti-diagonal vectorisation
        operates on cell slices, not per-cell visits.
        """
        if (ordinate_idx is None) == (direction_sign is None):
            raise ValueError(
                "dag_walk requires exactly one of `ordinate_idx` or "
                "`direction_sign`."
            )
        if self.reduced is None:
            raise ValueError(
                "dag_walk is only defined for meshes with a "
                "ReducedStreamingOperator (1-D Cartesian, spherical, "
                "or cylindrical).  2-D Cartesian wavefront sweeps "
                "use anti-diagonal scheduling, not per-cell visits."
            )
        coord = self.reduced.coord

        # Direction-keyed branch: resolve a non-degenerate representative
        # ordinate, then delegate to the ordinate-keyed branch (single source
        # of truth — Pattern 2). Cylindrical ``mu_level_idx`` is required by
        # both ``_representative_ordinate`` and the cylindrical visit iterator,
        # each via ``_require_mu_level`` (fail-loud at point of use).
        if direction_sign is not None:
            if direction_sign not in (+1, -1):
                raise ValueError(
                    f"direction_sign must be +1 or -1; got "
                    f"{direction_sign}"
                )
            ordinate_idx = self._representative_ordinate(
                direction_sign, mu_level_idx,
            )
        # Exactly one of ordinate_idx / direction_sign was supplied (the XOR
        # guard above); the direction-keyed arm resolved it, so the ordinate
        # is now concrete — narrow for the type checker.
        if ordinate_idx is None:  # pragma: no cover — unreachable per XOR guard
            raise ValueError(
                "dag_walk: ordinate_idx unresolved after mode dispatch."
            )

        # Ordinate-keyed branch.
        if coord is CoordSystem.CARTESIAN:
            yield from self._iter_cartesian_visits(ordinate_idx)
            return
        if coord is CoordSystem.SPHERICAL:
            yield from self._iter_spherical_visits(ordinate_idx)
            return
        if coord is CoordSystem.CYLINDRICAL:
            yield from self._iter_cylindrical_visits(
                ordinate_idx, self._require_mu_level(mu_level_idx),
            )
            return
        raise ValueError(  # pragma: no cover — exhaustive match above
            f"Unknown coord system: {coord!r}"
        )

    def dag_walk_cell_indices(
        self,
        *,
        direction_sign: int,
        mu_level_idx: int | None = None,
    ) -> Iterator[int]:
        r"""Lightweight twin of :meth:`dag_walk` — yields just cell indices.

        Consumers that build their own per-cell algebra from primitives
        (the loss-representation walk in
        :mod:`orpheus.sn.loss_representation` — the former unified matvec
        ``transport_operator_matvec_unified`` was its Depth-B predecessor,
        deleted at the walk unification) only need the cell traversal
        order, not the full
        :class:`~orpheus.transport.spatial.scheme.CellVisit` packet.

        Eliminates per-cell-per-call ``ReducedStreamingOperator.streaming_terms()``
        construction + frozen-dataclass overhead.  PR-TYPED-6c profiling
        showed this was ~14% of matvec time on slab, ~18% on cylinder
        — all building a packet the matvec discards.

        Cell-iteration order matches :meth:`dag_walk`:

        * Slab, sphere, cylinder non-degenerate: ``range(nx)`` for
          :math:`\mu_n \ge 0`, ``range(nx-1, -1, -1)`` for :math:`\mu_n < 0`.
        * Cylindrical pure-azimuthal degenerate
          (:math:`|\eta_n| < 10^{-15}`): ``range(nx)`` regardless of
          ``direction_sign`` — same as :meth:`dag_walk`.
        """
        if self.reduced is None:
            raise ValueError(
                "dag_walk_cell_indices is only defined for meshes with a "
                "ReducedStreamingOperator (1-D Cartesian, spherical, "
                "or cylindrical)."
            )
        coord = self.reduced.coord
        if direction_sign not in (+1, -1):
            raise ValueError(
                f"direction_sign must be +1 or -1; got {direction_sign}"
            )

        # Resolve the representative ordinate's signed primary cosine.
        # (Cylindrical mu_level_idx is required by _representative_ordinate
        # and the global-ordinate lookup below, each via _require_mu_level.)
        ordinate_idx = self._representative_ordinate(
            direction_sign, mu_level_idx,
        )
        if coord is CoordSystem.CYLINDRICAL:
            mu_level = self._require_mu_level(mu_level_idx)
            level_indices = self.quad.level_indices  # type: ignore[attr-defined]
            global_n = int(level_indices[mu_level][ordinate_idx])
            mu_n = float(self.quad.mu_x[global_n])
        else:
            mu_n = float(self.quad.mu_x[ordinate_idx])

        # Cylindrical degenerate ordinates iterate forward regardless of sign.
        if (
            coord is CoordSystem.CYLINDRICAL
            and abs(mu_n) < self._DEGENERATE_ABS_ETA_THRESHOLD
        ):
            yield from range(self.nx)
            return

        if mu_n >= 0:
            yield from range(self.nx)
        else:
            yield from range(self.nx - 1, -1, -1)

    def _require_mu_level(self, mu_level_idx: int | None) -> int:
        """Narrow ``mu_level_idx`` to ``int`` for a cylindrical sweep.

        Cylindrical 1-D radial sweeps are organised by μ-level (a subset of
        azimuthal ordinates at one polar cosine), so every cylindrical
        traversal needs ``mu_level_idx``; slab/sphere pass ``None``. Single
        source of truth for the "cylindrical requires mu_level_idx" contract
        (Pattern 2) — fails loudly (``-O``-safe, not ``assert``) and returns
        the narrowed ``int`` so callers index ``level_indices`` cleanly.
        """
        if mu_level_idx is None:
            raise ValueError(
                "cylindrical sweep requires mu_level_idx (which μ-level the "
                "ordinate subset belongs to); slab/sphere pass None."
            )
        return mu_level_idx

    def _representative_ordinate(
        self,
        direction_sign: int,
        mu_level_idx: int | None,
    ) -> int:
        """Pick a non-degenerate ordinate matching the direction sign.

        Cell ordering in :meth:`dag_walk` depends only on
        ``direction_sign`` (and the level for cylindrical), so any
        non-degenerate ordinate in the correct sign class produces
        the same cell sequence.  The degenerate :math:`|\\eta| <
        10^{-15}` ordinates are excluded because they iterate forward
        regardless of sign and would not match the bulk direction's
        signed iteration.
        """
        assert self.reduced is not None
        coord = self.reduced.coord
        eps = self._DEGENERATE_ABS_ETA_THRESHOLD
        if coord is CoordSystem.CYLINDRICAL:
            mu_level = self._require_mu_level(mu_level_idx)
            level_indices = self.quad.level_indices  # type: ignore[attr-defined]
            level_ords = np.asarray(level_indices[mu_level])
            eta_at_level = self.quad.eta[level_ords]
            if direction_sign == +1:
                cand = np.where(eta_at_level > +eps)[0]
            else:
                cand = np.where(eta_at_level < -eps)[0]
            if cand.size == 0:
                raise ValueError(
                    f"No non-degenerate ordinate in cylindrical level "
                    f"{mu_level} satisfies "
                    f"direction_sign={direction_sign}."
                )
            return int(cand[0])
        mu_x = self.quad.mu_x
        if direction_sign == +1:
            cand = np.where(mu_x > +eps)[0]
        else:
            cand = np.where(mu_x < -eps)[0]
        if cand.size == 0:
            raise ValueError(
                f"No non-degenerate ordinate satisfies "
                f"direction_sign={direction_sign} in this quadrature."
            )
        return int(cand[0])

    def _make_cell_visit(
        self,
        *,
        cell_idx: int,
        global_ordinate: int,
        face_area_downstream: float,
        st: StreamingTerms,
    ) -> CellVisit:
        r"""Assemble one :class:`CellVisit`, sourcing the angular-closure τ / c.

        Issue #236 Phase 2 B2/B3 — the single production site that stamps the
        Morel--Montry angular weight :attr:`CellVisit.tau` (B3) and the
        derived weighted-diamond constants
        (:attr:`CellVisit.c_in` / :attr:`CellVisit.c_out`, B2) onto a visit.
        ALL four ``dag_walk`` yield paths (slab / sphere / cylinder /
        cylindrical-degenerate) funnel through here so the lookup lives
        in exactly ONE place (Pattern 2 — no per-site divergence).

        The values are read from the mesh's canonical angular-closure
        owner :attr:`pole_angular_closure` via its per-global-ordinate
        ``(N,)`` accessors
        (:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
        / ``c_in_per_ordinate`` / ``c_out_per_ordinate``) — NOT rebuilt from
        ``st.alpha_*`` / ``st.tau_mm`` (the inline formula the former
        duplication sites carried).  ``global_ordinate`` is the
        GLOBAL ordinate index: ``direction_idx`` for slab / sphere,
        ``level_indices[mu_level_idx][m]`` for cylinder (mirroring the
        index :meth:`streaming_terms` resolves).  Slab / Cartesian reads
        the identity closure's neutral values, so ``c_in == c_out == 0.0``
        and ``tau == 1.0`` there.
        """
        closure = self.pole_angular_closure
        return CellVisit(
            cell_idx=cell_idx,
            streaming_terms=st,
            face_area_downstream=face_area_downstream,
            c_in=float(closure.c_in_per_ordinate[global_ordinate]),
            c_out=float(closure.c_out_per_ordinate[global_ordinate]),
            tau=float(closure.tau_per_ordinate[global_ordinate]),
        )

    def _iter_cartesian_visits(
        self,
        ordinate_idx: int,
    ) -> Iterator[CellVisit]:
        """Yield slab (1-D Cartesian) visits in sweep direction.

        Order: forward (cell 0 → nx-1) for :math:`\\mu \\ge 0`,
        backward for :math:`\\mu < 0`.  Slab carries
        ``face_area_downstream = 1.0`` (neutral curvature; Issue
        #196 Phase G Step 2.5) so the unified cell-balance helper
        consumes one geometry-blind number.
        """
        assert self.reduced is not None
        mu_n = float(self.quad.mu_x[ordinate_idx])
        cell_indices = (
            range(self.nx) if mu_n >= 0 else range(self.nx - 1, -1, -1)
        )
        for i in cell_indices:
            st = self.reduced.streaming_terms(
                cell_idx=i, direction_idx=ordinate_idx,
            )
            # Slab: ``direction_idx`` IS the global ordinate; the identity
            # closure yields neutral c == 0.0 (#236 Phase 2 B2).
            yield self._make_cell_visit(
                cell_idx=i,
                global_ordinate=ordinate_idx,
                face_area_downstream=1.0,
                st=st,
            )

    def _iter_spherical_visits(
        self,
        ordinate_idx: int,
    ) -> Iterator[CellVisit]:
        """Yield spherical visits in sweep direction.

        Outward (:math:`\\mu \\ge 0`): cell 0 → nx-1, downstream face
        is the outer face ``A[i+1]``.  Inward (:math:`\\mu < 0`):
        cell nx-1 → 0, downstream face is the inner face ``A[i]``.
        """
        assert self.reduced is not None
        mu_n = float(self.quad.mu_x[ordinate_idx])
        if mu_n >= 0:
            cell_indices = range(self.nx)
            select_outer = True
        else:
            cell_indices = range(self.nx - 1, -1, -1)
            select_outer = False
        for i in cell_indices:
            st = self.reduced.streaming_terms(
                cell_idx=i, direction_idx=ordinate_idx,
            )
            face_downstream = (
                st.face_area_outer if select_outer else st.face_area_inner
            )
            # Sphere: ``direction_idx`` IS the global ordinate (#236 B2).
            yield self._make_cell_visit(
                cell_idx=i,
                global_ordinate=ordinate_idx,
                face_area_downstream=face_downstream,
                st=st,
            )

    def _iter_cylindrical_visits(
        self,
        ordinate_idx: int,
        mu_level_idx: int,
    ) -> Iterator[CellVisit]:
        """Yield cylindrical visits in sweep direction for one level.

        ``ordinate_idx`` is the within-level azimuthal index
        :math:`m \\in [0, M)`.  The global ordinate is resolved via
        ``quad.level_indices[mu_level_idx][ordinate_idx]``.

        * :math:`\\eta_n \\ge 0` outward: cell 0 → nx-1, downstream
          face is the outer face.
        * :math:`\\eta_n < 0` inward: cell nx-1 → 0, downstream
          face is the inner face.
        * :math:`|\\eta_n| < 10^{-15}` pure-azimuthal degenerate:
          forward iteration (so the angular M-M closure runs in a
          natural order) but ``face_area_downstream`` is ``None`` —
          no spatial face flow.
        """
        assert self.reduced is not None
        level_indices = self.quad.level_indices  # type: ignore[attr-defined]
        global_n = int(level_indices[mu_level_idx][ordinate_idx])
        eta_n = float(self.quad.eta[global_n])
        abs_eta = abs(eta_n)

        if abs_eta < self._DEGENERATE_ABS_ETA_THRESHOLD:
            # Pure-azimuthal degenerate: no spatial flow.  Iterate
            # forward so the angular M-M closure runs in a natural
            # order; ``face_area_downstream = 0.0`` signals "no
            # spatial flow" to the strategy (geometric truth — the
            # cell has no radial face on this ordinate).  Issue #196
            # Phase G Step 2.5: replaced ``None`` with the
            # geometrically-correct float ``0.0``.
            for i in range(self.nx):
                st = self.reduced.streaming_terms(
                    cell_idx=i,
                    direction_idx=ordinate_idx,
                    mu_level_idx=mu_level_idx,
                )
                # Cylinder: the global ordinate is resolved through the
                # level partition (``global_n`` above) — the SAME index
                # ``streaming_terms`` used to read the per-level α / τ
                # (#236 Phase 2 B2).
                yield self._make_cell_visit(
                    cell_idx=i,
                    global_ordinate=global_n,
                    face_area_downstream=0.0,
                    st=st,
                )
            return

        if eta_n >= 0:
            cell_indices = range(self.nx)
            select_outer = True
        else:
            cell_indices = range(self.nx - 1, -1, -1)
            select_outer = False
        for i in cell_indices:
            st = self.reduced.streaming_terms(
                cell_idx=i,
                direction_idx=ordinate_idx,
                mu_level_idx=mu_level_idx,
            )
            face_downstream = (
                st.face_area_outer if select_outer else st.face_area_inner
            )
            # Cylinder: global ordinate via the level partition (#236 B2).
            yield self._make_cell_visit(
                cell_idx=i,
                global_ordinate=global_n,
                face_area_downstream=face_downstream,
                st=st,
            )

    # ── Stencil setup ─────────────────────────────────────────────────

    def _setup_cartesian(self) -> None:
        r"""Precompute the raw per-axis down-face streaming coefficient ``g``.

        The **scheme-agnostic** geometric streaming, one array per spatial
        axis:

        .. math::

            g_a \;=\; \frac{|\mu_a|\,A_{\rm down}}{V}
                \;=\; \frac{|\mu_a|}{\Delta a}
                \qquad(\text{Cartesian tensor-product: } A_{\rm down}/V = 1/\Delta a)

        This is NOT the DD denominator term.  The diamond-difference closure
        contributes :math:`\Sigma_t + \sum_a 2g_a` to the cell-balance
        denominator — the factor :math:`2 = 1/w_{\rm DD}` is DD's diamond
        closure (:math:`\psi_{\rm out} = 2\bar\psi - \psi_{\rm in}`), owned by
        the *scheme* and applied inside its cell kernel, NOT baked into this
        geometric accessor (#240).  Linear-Discontinuous reads the same raw
        ``g`` without DD's factor.

        Precomputing avoids per-ordinate per-cell divisions in the inner
        sweep loop.  Built over ``range(ndim)`` from the canonical per-axis
        accessors (``quad.axis_cosines(a)`` — the legacy ``mu_x`` / ``mu_y``
        names are property views of exactly these columns), with NO phantom
        axis: a 1-D mesh carries one streaming array, not an ``ny=1`` second.
        """
        # _streaming_axes[a][n, i] = |μ_a[n]| / Δa[i] — the RAW down-face
        # streaming g (shape (N_ord, n_a)); the scheme owns its closure factor.
        self._streaming_axes: tuple[np.ndarray, ...] | None = tuple(
            np.abs(self.quad.axis_cosines(a))[:, None]
            / widths[None, :]
            for a, widths in enumerate(self.axis_widths)
        )

        # Curvature terms (None for Cartesian — placeholder for curvilinear)
        self.curvature = None

    # ── Backward-compat property accessors ────────────────────────────
    #
    # These properties route to ``self.reduced`` (the
    # :class:`ReducedStreamingOperator` built at construction) and emit a
    # ``DeprecationWarning`` — use ``self.reduced.<name>`` directly.

    @property
    def face_areas(self) -> np.ndarray:
        """[Deprecated] Cell face areas. Use ``self.reduced.face_areas`` instead."""
        warnings.warn(
            "SNMesh.face_areas is deprecated; "
            "use SNMesh.reduced.face_areas instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.face_areas is not None
        return self.reduced.face_areas

    @property
    def delta_A(self) -> np.ndarray:
        """[Deprecated] Face-area differences. Use ``self.reduced.delta_A`` instead."""
        warnings.warn(
            "SNMesh.delta_A is deprecated; "
            "use SNMesh.reduced.delta_A instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        assert self.reduced.delta_A is not None
        return self.reduced.delta_A

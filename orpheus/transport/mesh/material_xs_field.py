r"""Macroscopic cross-section field over an SN domain.

Issue #197 PR-TYPED-1 — the typed wrapper that closes the 8 leaked
per-material dispatch loops scattered across :mod:`orpheus.sn.scattering`
and :mod:`orpheus.sn.solver`.  Before this PR the per-material structure
of the cross-section data leaked through every consumer:

* :class:`~orpheus.sn.scattering.ScatteringOperator` carried a
  ``cells_by_mat: dict[int, (ix, iy)]`` constructor parameter and
  iterated it explicitly in :meth:`~ScatteringOperator.add_iso_source`,
  :meth:`~ScatteringOperator.add_n2n_source`,
  :meth:`~ScatteringOperator.foldable_part`,
  :meth:`~ScatteringOperator.residual_part`,
  :meth:`~ScatteringOperator.is_foldable_into_sigma_r`,
  :meth:`~ScatteringOperator.foldable_sigma`, and
  :class:`~orpheus.sn.scattering.LegendreMomentScattering.apply`.
* :class:`~orpheus.sn.solver.SNSolver` carried a parallel
  ``_cells_by_mat`` and seven separate XS attributes
  (``sig_t``, ``sig_a``, ``sig_p``, ``chi``, ``sig_s``, ``sig2``,
  ``sig_s0``) all keyed on the same per-material/per-cell topology.
* :meth:`~SNSolver.compute_group_production_rate` ran an explicit
  ``for mid, (ix, iy) in self._cells_by_mat.items()`` to assemble
  the (n,2n) contribution.

The single source of truth is :class:`MaterialXSField`, a typed wrapper
over the per-material :class:`~orpheus.data.macro_xs.mixture.Mixture`
dict plus the spatial distribution carried by the :class:`SNMesh`'s
``mat_map``.  Two access modes:

* **Per-cell views** (``total_cross_section``, ``absorption_cross_section``,
  ``fission_production``, ``emission_spectrum``) — the broadcast
  ``(ng, nx, ny)`` arrays every operator's per-cell math consumes.
  Built lazily on first access via :func:`assemble_cell_xs`; cached.
* **Per-material accessors** (``scattering_legendre``, ``n2n_matrix``,
  ``fission_production_per_material``, ``chi_per_material``,
  ``cells_by_material``) — for operations that genuinely exploit
  per-material structure (group-coupling matmul on small ``(ng, ng)``
  matrices, anisotropic moment scattering).
* **Named typed verbs** (``apply_p0_in_scatter``, ``apply_n2n``,
  ``apply_legendre_scattering_moments``, ``add_n2n_to_group_rate``)
  — the lifted-up forms of the formerly-leaked per-material dispatch
  loops.  Each verb captures one piece of math and reads as the
  domain (``coding-elegance`` Pattern 1 + Pattern 3).

Composability framing (the user's three operations):

* **Mixing**: weighted volume-average of two ``MaterialXSField``\ s
  → single homogenised ``MaterialXSField``.  Not yet implemented;
  the homogenisation step lives outside SN today.  Future Wave.
* **Restriction**: per-region subset returns a ``MaterialXSField``
  on the smaller domain.  Future Wave (CP / MoC consumer pattern).
* **Action**: ``mat_xs.fission_production * scalar_flux`` eventually
  reads as the math (after PR-TYPED-2 introduces typed flux fields
  with ``__mul__`` dunders).  This PR keeps bare ``np.ndarray``
  inputs/outputs; the typed-action wiring lands in PR-TYPED-2.

Storage discipline (``coding-elegance`` Pattern 7 — normalise at
definition site): the per-material :class:`Mixture` dict + the
:class:`SNMesh` ARE the source of truth.  The lazy per-cell views
are cached but content-derived from the frozen inputs; they cannot
diverge from the source data.

Capability matrix (which sites collapsed):

==========  ============================================  =================================
Old site    Old pattern                                   New call
==========  ============================================  =================================
scattering  ``for mid in cells_by_mat: ... add ...``      :meth:`apply_p0_in_scatter`
scattering  ``for mid in cells_by_mat: ... 2 · sig2 ...`` :meth:`apply_n2n`
scattering  ``for mid in cells_by_mat: sig_s[mid][l] @ M``:meth:`apply_legendre_scattering_moments`
scattering  ``for mid in sig_s.items(): np.diag(...)``    :meth:`foldable_sig_s` (helper)
scattering  ``for mid in sig_s.items(): off-diag``        :meth:`residual_sig_s` (helper)
scattering  ``for mid in sig_s.items(): np.allclose(...)``:meth:`is_p0_diagonal_with_zero_n2n`
scattering  ``for mid in sig_s.items(): diag(...)``       :meth:`foldable_sigma`
solver      ``for mid in _cells_by_mat: sig2 ...``        :meth:`add_n2n_to_group_rate`
==========  ============================================  =================================

Units (the discipline that physics code make units explicit, per
``coding-elegance`` Pattern 3): macroscopic cross sections in
``1/cm`` per energy group; ``chi`` is dimensionless emission spectrum;
the per-material/per-cell expansion is broadcast across cells, units
unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from orpheus.data.macro_xs.cell_xs import assemble_cell_xs
from orpheus.transport.fields.cross_section_field import CrossSectionField

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.numerics.frame import FrameBase
    # The ``mesh`` field is typed against ``SNMesh`` under TYPE_CHECKING
    # only (an L2-field / L3-mesh annotation, layer-legal — the runtime
    # never imports sn). MaterialXSField in fact reads only
    # ``MaterialMesh`` data (``materials`` / ``mat_map`` / ``ng`` /
    # ``spatial_shape``); retyping the whole typed-field hierarchy's
    # ``mesh`` annotation from ``SNMesh`` to ``MaterialMesh`` (splitting
    # the bulk-data ``ScalarField`` branch from the quad/trace-dependent
    # ``AngularField`` / ``BoundaryField`` branches) is a separate L2
    # follow-up — see issue #267.
    from orpheus.sn.geometry import SNMesh


__all__ = ["MaterialXSField"]


@dataclass
class MaterialXSField:
    r"""Macroscopic cross-section field over an SN phase space.

    Owns the per-material :class:`~orpheus.data.macro_xs.mixture.Mixture`
    data plus the spatial distribution via the SN mesh's ``mat_map``.
    Exposes BOTH per-material accessors (for operations that exploit
    per-material structure, e.g. group-coupling matmul on ``(ng, ng)``
    matrices) AND per-cell expanded views (for operations that need
    cell-grid layout, e.g. the streaming/collision algebra).

    The per-cell views (``total_cross_section`` etc.) are CACHED on
    first access via :func:`~orpheus.data.macro_xs.cell_xs.assemble_cell_xs`;
    the per-material side carries the source of truth.

    Parameters
    ----------
    materials : dict[int, Mixture]
        Per-material macroscopic cross sections keyed by integer
        material id.  All materials must agree on ``ng`` (already
        validated by :class:`SNMesh` at construction).
    mesh : SNMesh
        The SN phase-space carrier — supplies ``mat_map``, ``nx``,
        ``ny``, ``ng``, plus the geometry/quadrature handles that
        downstream consumers reach through ``mat_xs.mesh.*``.

    Attributes (cached)
    -------------------
    All lazy: populated on first read of the corresponding property
    via :meth:`__post_init__`-free direct construction.  Mutating
    the underlying :attr:`materials` / :attr:`mesh` after construction
    is undefined behaviour — :class:`MaterialXSField` is conceptually
    frozen but uses a non-frozen dataclass to make the lazy caches
    natural.

    Notes
    -----
    Non-frozen dataclass (not :pyfunc:`dataclasses.dataclass(frozen=True)`)
    because the lazy per-cell caches mutate ``self``.  The
    frozen-with-``object.__setattr__`` workaround would obscure the
    storage discipline and complicate testing.  The "frozenness" we
    want is content-immutability of :attr:`materials` and :attr:`mesh`,
    which Python doesn't enforce structurally but which every consumer
    of this class respects by convention.
    """

    materials: dict[int, "Mixture"]
    mesh: "SNMesh"

    # Lazy per-cell views — populated on first access.
    _sig_t_cell: np.ndarray | None = field(default=None, init=False, repr=False)
    _sig_a_cell: np.ndarray | None = field(default=None, init=False, repr=False)
    _sig_p_cell: np.ndarray | None = field(default=None, init=False, repr=False)
    _chi_cell: np.ndarray | None = field(default=None, init=False, repr=False)
    _cells_by_mat: dict[int, tuple[np.ndarray, np.ndarray]] | None = field(
        default=None, init=False, repr=False,
    )
    # Cached dense (n,2n) matrices to avoid repeated ``.todense()``
    # in the (n,2n) hot path.
    _n2n_dense: dict[int, np.ndarray] | None = field(
        default=None, init=False, repr=False,
    )
    # Cached per-material Legendre scattering lists (already dense).
    _sig_s_dense: dict[int, list[np.ndarray]] | None = field(
        default=None, init=False, repr=False,
    )

    # ── Construction helpers ──────────────────────────────────────────

    @classmethod
    def from_mesh(cls, mesh: "SNMesh") -> "MaterialXSField":
        """Build the XS field directly from the mesh's authoritative materials.

        Standard constructor — the materials dict already lives on the
        mesh (Issue #197 PR-TYPED-0).  Reads only
        :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` data
        (``mesh.materials`` / ``mesh.mat_map`` / ``mesh.ng`` /
        ``mesh.spatial_shape``), so it accepts any ``MaterialMesh`` or
        subclass (``SNMesh``); the parameter is typed ``SNMesh`` only
        pending the typed-field ``mesh`` retype (see module note).
        """
        return cls(materials=mesh.materials, mesh=mesh)

    @classmethod
    def _synthetic_for_tests(
        cls,
        *,
        sig_s: dict[int, list[np.ndarray]],
        sig2: dict[int, np.ndarray] | None = None,
        cells_by_mat: dict[int, tuple[np.ndarray, np.ndarray]],
        ng: int,
        nx: int,
        ny: int = 1,
    ) -> "MaterialXSField":
        """Build a synthetic :class:`MaterialXSField` from raw dicts.

        Test-only construction path that bypasses :class:`SNMesh` and
        the :class:`Mixture` registry.  Lets foundation tests for
        :class:`LegendreMomentScattering` and
        :class:`ScatteringOperator` (which exercise the per-material
        dispatch in isolation) build the typed wrapper without paying
        the full mesh-construction cost.

        The returned field carries pre-populated dense caches; it has
        NO :attr:`mesh` reference (set to a minimal stand-in carrying
        ``nx``, ``ny``, ``ng``, and ``mat_map`` derived from
        ``cells_by_mat``).  Consumers that access per-cell views
        (:attr:`total_cross_section` etc.) will fail — those views
        require real per-material :class:`Mixture` data.  This factory
        is for the per-material-dispatch path only.

        Parameters
        ----------
        sig_s : dict[int, list[np.ndarray]]
            Per-material Legendre scattering lists.  Each entry is a
            list of ``(ng, ng)`` arrays.
        sig2 : dict[int, np.ndarray] or None
            Per-material ``(n,2n)`` matrices.  Defaults to zero when
            omitted.
        cells_by_mat : dict[int, tuple[np.ndarray, np.ndarray]]
            Per-material ``(ix, iy)`` cell-index arrays.
        ng, nx, ny : int
            Phase-space sizes.

        Returns
        -------
        MaterialXSField
            Synthetic instance suitable for foundation tests of the
            per-material dispatch code paths.
        """
        # Minimal mesh stand-in carrying only the attributes the
        # per-material verbs reach for.
        mat_map = np.full((nx, ny), -1, dtype=int)
        for mid, (ix, iy) in cells_by_mat.items():
            mat_map[ix, iy] = mid

        class _StubMesh:
            """Tiny stand-in covering only what MaterialXSField reads."""
            def __init__(self) -> None:
                self.spatial_shape = (nx, ny)
                self.ng = ng
                self.mat_map = mat_map
                # The materials dict — minimal Mixture-shaped stand-in.
                # MaterialXSField.materials reads ``.SigS``, ``.Sig2``,
                # ``.SigP``, ``.chi`` from each entry; for the
                # per-material dispatch verbs only ``sig_s`` and
                # ``sig2`` matter.

        # Materials dict: only the keys matter (every per-material
        # accessor routes through _sig_s_dense / _n2n_dense, which we
        # pre-populate below).
        materials: dict[int, "Mixture"] = {
            mid: None  # type: ignore[dict-item]
            for mid in sig_s
        }
        # Build dense caches from inputs.
        sig_s_dense = {
            mid: [np.asarray(s) for s in mats]
            for mid, mats in sig_s.items()
        }
        if sig2 is None:
            n2n_dense = {
                mid: np.zeros((ng, ng)) for mid in sig_s
            }
        else:
            n2n_dense = {mid: np.asarray(m) for mid, m in sig2.items()}

        out = cls(materials=materials, mesh=_StubMesh())  # type: ignore[arg-type]
        out._sig_s_dense = sig_s_dense
        out._n2n_dense = n2n_dense
        out._cells_by_mat = cells_by_mat
        return out

    def with_overridden_sig_s_and_n2n(
        self,
        sig_s_dense: dict[int, list[np.ndarray]],
        n2n_dense: dict[int, np.ndarray],
    ) -> "MaterialXSField":
        r"""Return a sibling carrying overridden scattering / (n,2n) data.

        Used by :meth:`ScatteringOperator.foldable_part` and
        :meth:`ScatteringOperator.residual_part` to build derived
        scattering operators where the per-material P0 / Pℓ / (n,2n)
        matrices are modified (diagonal-only / off-diagonal split)
        without rebuilding the full :class:`Mixture` dict.

        The sibling shares :attr:`mesh` and the cell-grid views
        (:attr:`total_cross_section`, ``etc.``) with ``self`` — the
        derived scattering data only affects per-material accessors
        and the typed verbs that consume them.

        Parameters
        ----------
        sig_s_dense : dict[int, list[np.ndarray]]
            New per-material Legendre scattering list (one ``(ng, ng)``
            matrix per Legendre order).  Keys must match
            :attr:`materials`.
        n2n_dense : dict[int, np.ndarray]
            New per-material ``(n,2n)`` matrix.  Keys must match
            :attr:`materials`.

        Returns
        -------
        MaterialXSField
            Sibling with overridden caches populated.
        """
        sibling = MaterialXSField(materials=self.materials, mesh=self.mesh)
        # Pre-populate the dense caches with the overrides; the
        # per-cell views and cells_by_material will be lazily
        # populated on first access (identical to self's by design).
        sibling._sig_s_dense = sig_s_dense
        sibling._n2n_dense = n2n_dense
        return sibling

    # ── Homogenisation: project the whole field through coarse frames ──

    def project_through(
        self, sigma_frame: "FrameBase", emission_frame: "FrameBase", /,
    ) -> dict[int, "Mixture"]:
        r"""Homogenise the whole cross-section field through two coarse frames.

        Collapse every per-fine-cell channel onto the coarse cells of the frames'
        shared trial (cell-indicator) basis, returning one effective
        :class:`~orpheus.data.macro_xs.mixture.Mixture` per coarse cell. The field
        owns the **channel → weighting taxonomy** and routes accordingly:

        * the **rate-bearing** channels — :math:`\Sigma_t,\Sigma_c,\Sigma_L,
          \Sigma_f,\nu\Sigma_f` (vectors) and :math:`\Sigma_{s,\ell},\Sigma_{2n}`
          (``[g_from, g_to]`` matrices) — collapse through ``sigma_frame``, whose
          flux-weighted test basis makes :meth:`~orpheus.numerics.frame.FrameBase.project`
          the reaction-rate-preserving average :math:`\Sigma_R = \int_R\varphi\Sigma\,
          \mathrm{d}V/\int_R\varphi\,\mathrm{d}V` (matrices weight by the **source**
          group, the leading axis the test weight aligns to);
        * the **emission spectrum** :math:`\chi` collapses through ``emission_frame``,
          whose production-weighted test (:math:`p=\sum_g\nu\Sigma_{f,g}\varphi_g`)
          gives the production-weighted convex average :math:`\chi_R = \int_R p\chi\,
          \mathrm{d}V/\int_R p\,\mathrm{d}V` (a convex combination of simplices, hence
          a simplex).

        Two frames because the two collapses preserve two *different* conserved
        functionals (reaction rate vs emission rate) — the campaign's "one frame
        carries one test weighting" ruling. Both frames share the trial basis +
        geometric measure; the caller (which owns the flux) builds their test
        weightings.

        Parameters
        ----------
        sigma_frame : FrameBase
            The flux-weighted homogenisation frame for the rate-bearing channels.
        emission_frame : FrameBase
            The production-weighted frame for :math:`\chi`.

        Returns
        -------
        dict[int, Mixture]
            One effective :class:`Mixture` per coarse cell, keyed by coarse-cell
            index ``0 .. n_coarse-1``. A group with zero region weight collapses to a
            zero effective cross section there (the frame's Moore–Penrose Gram
            pseudo-inverse — no reaction rate to preserve).
        """
        from scipy.sparse import csr_matrix

        from orpheus.data.macro_xs.mixture import Mixture

        materials = self.materials
        ng = self.ng
        mat_of_fine = np.asarray(self.mesh.mat_map, dtype=int).ravel()

        def _gather_vector(attr: str) -> np.ndarray:
            """Per-fine-cell view of a ``(ng,)`` Mixture channel — ``(n_fine, ng)``."""
            return np.array([getattr(materials[m], attr) for m in mat_of_fine])

        n_legendre = max(len(materials[m].SigS) for m in materials)

        def _gather_legendre(order: int) -> np.ndarray:
            """Per-fine-cell dense ``Σ_{s,ℓ}`` — ``(n_fine, ng, ng)``; zero-pad short lists."""
            return np.array([
                np.asarray(materials[m].SigS[order].todense())
                if order < len(materials[m].SigS) else np.zeros((ng, ng))
                for m in mat_of_fine
            ])

        # Rate-bearing channels → the flux-weighted frame. ``project`` IS the
        # rate-preserving collapse G⁻¹M (the per-channel inline gather/collapse the
        # method body used to carry now lives here, projected as ONE field).
        sig_t = sigma_frame.project(_gather_vector("SigT"))
        sig_c = sigma_frame.project(_gather_vector("SigC"))
        sig_l = sigma_frame.project(_gather_vector("SigL"))
        sig_f = sigma_frame.project(_gather_vector("SigF"))
        sig_p = sigma_frame.project(_gather_vector("SigP"))
        sig_s = [sigma_frame.project(_gather_legendre(l)) for l in range(n_legendre)]
        sig2 = sigma_frame.project(
            np.array([np.asarray(materials[m].Sig2.todense()) for m in mat_of_fine])
        )

        # χ → the production-weighted frame (a different conserved rate).
        chi = emission_frame.project(_gather_vector("chi"))

        eg = next(iter(materials.values())).eg
        n_coarse = sig_t.shape[0]
        return {
            region: Mixture(
                SigC=sig_c[region], SigL=sig_l[region], SigF=sig_f[region],
                SigP=sig_p[region], SigT=sig_t[region],
                SigS=[csr_matrix(sig_s[l][region]) for l in range(n_legendre)],
                Sig2=csr_matrix(sig2[region]), chi=chi[region], eg=eg,
            )
            for region in range(n_coarse)
        }

    # ── Per-cell views (lazy, cached) ─────────────────────────────────

    @property
    def total_cross_section(self) -> np.ndarray:
        r""":math:`\sigma_t` per-cell view, shape ``(ng, nx, ny)``.

        Units: ``1/cm`` per energy group.  Built lazily on first
        access via :func:`~orpheus.data.macro_xs.cell_xs.assemble_cell_xs`
        with a ``.T.reshape(ng, nx, ny)`` to the principled layout
        (Issue #196 PR-INDEX-3).  Cached.
        """
        if self._sig_t_cell is None:
            self._ensure_cell_views()
        return self._sig_t_cell  # type: ignore[return-value]

    @property
    def absorption_cross_section(self) -> np.ndarray:
        r""":math:`\sigma_a` per-cell view, shape ``(ng, nx, ny)``."""
        if self._sig_a_cell is None:
            self._ensure_cell_views()
        return self._sig_a_cell  # type: ignore[return-value]

    @property
    def fission_production(self) -> np.ndarray:
        r""":math:`\nu \Sigma_f` per-cell view, shape ``(ng, nx, ny)``.

        Production cross-section ``νΣ_f`` — the rate at which fission
        emits new neutrons per absorption per group.  Units: ``1/cm``.
        """
        if self._sig_p_cell is None:
            self._ensure_cell_views()
        return self._sig_p_cell  # type: ignore[return-value]

    @property
    def emission_spectrum(self) -> np.ndarray:
        r""":math:`\chi` per-cell view, shape ``(ng, nx, ny)``.

        Fission emission spectrum (dimensionless; ``Σ_g χ_g = 1`` per
        fissile material).  Broadcast across cells.
        """
        if self._chi_cell is None:
            self._ensure_cell_views()
        return self._chi_cell  # type: ignore[return-value]

    def _ensure_cell_views(self) -> None:
        """Populate the four per-cell views via :func:`assemble_cell_xs`.

        Single producer of the principled-layout per-cell arrays.
        Run-once: subsequent accesses hit the cache.  This is the
        canonical bridge between the per-material :class:`Mixture`
        flat representation and the principled per-cell ``(ng, nx, ny)``
        layout.
        """
        xs = assemble_cell_xs(self.materials, self.mesh.mat_map)
        ng = self.mesh.ng
        spatial = self.mesh.spatial_shape
        # .T.reshape: producer emits (N_cells, ng); flip to (ng, N_cells)
        # then split N_cells back into (*spatial) — the principled
        # (ng, *spatial) layout (rank == ndim; no phantom ny=1 on 1-D).
        self._sig_t_cell = xs.sig_t.T.reshape(ng, *spatial)
        self._sig_a_cell = xs.sig_a.T.reshape(ng, *spatial)
        self._sig_p_cell = xs.sig_p.T.reshape(ng, *spatial)
        self._chi_cell = xs.chi.T.reshape(ng, *spatial)

    # ── Typed per-cell views (the CrossSectionField promotion, #257 S2) ──
    #
    # The field side of the operator promotion C = M[σ_t] (#257 §5.7): each
    # macroscopic cross section, wrapped as the typed
    # :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
    # (units 1/cm; physical cross sections, so cone-valued — but the cone is a
    # property, NOT a per-field invariant: a signed σ−σ′ is still a valid
    # CrossSectionField). These wrap the SAME cached ndarray as the raw views
    # above (``.values is`` the raw array — bit-identical, zero copy), so they
    # are a pure-additive typed lens, NOT a second representation. The raw
    # ndarray views remain the live consumer path; the typed accessors are the
    # migration target the S3 ``MultiplicationOperator`` reads from.
    #
    # All three σ are typed for SYMMETRY even though only σ_t (S3 collision) and
    # νΣ_f (S6 fission) have near-term consumers: the wrap is one proven pattern
    # (its rule-of-two is met in-diff), and a two-typed-one-raw surface would
    # invite a future σ_a consumer to open-code the wrap inline (a single-source
    # violation). Do NOT delete ``absorption_cross_section_field`` as "unused".

    @property
    def total_cross_section_field(self) -> CrossSectionField:
        r""":math:`\Sigma_t` as a typed :class:`CrossSectionField` (1/cm)."""
        return CrossSectionField.from_mesh(self.total_cross_section, self.mesh)

    @property
    def absorption_cross_section_field(self) -> CrossSectionField:
        r""":math:`\Sigma_a` as a typed :class:`CrossSectionField` (1/cm)."""
        return CrossSectionField.from_mesh(self.absorption_cross_section, self.mesh)

    @property
    def fission_production_field(self) -> CrossSectionField:
        r""":math:`\nu\Sigma_f` as a typed :class:`CrossSectionField` (1/cm)."""
        return CrossSectionField.from_mesh(self.fission_production, self.mesh)

    # ── Per-material accessors (source of truth) ─────────────────────

    @property
    def cells_by_material(self) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        r"""Per-material ``(ix, iy)`` index arrays — cached.

        For each material id, returns the ``(ix_array, iy_array)``
        pair such that ``mat_map[ix, iy] == mid`` everywhere.  This
        is the single index map the formerly-leaked per-material
        dispatch loops keyed on.  Most consumers should NOT use this
        directly — call one of the typed verbs
        (:meth:`apply_p0_in_scatter`, :meth:`apply_n2n`, etc.) that
        encapsulates the loop.

        Returns
        -------
        dict[int, tuple[np.ndarray, np.ndarray]]
            ``{mid: (ix, iy)}`` for every material id in
            :attr:`materials`.
        """
        if self._cells_by_mat is None:
            self._cells_by_mat = {
                mid: np.where(self.mesh.mat_map == mid)
                for mid in self.materials
            }
        return self._cells_by_mat

    def sig_s_legendre(self, material_id: int) -> list[np.ndarray]:
        r"""Per-material list of dense Legendre scattering matrices.

        Returns ``[Σ_{s,0}, Σ_{s,1}, ..., Σ_{s,L}]`` for the requested
        material, each entry a dense ``(ng, ng)`` array indexed
        ``[g_from, g_to]``.  Materially equivalent to
        ``[mix.SigS[l].todense() for l in ...]`` but cached.
        """
        if self._sig_s_dense is None:
            self._build_dense_caches()
        return self._sig_s_dense[material_id]  # type: ignore[index]

    def n2n_matrix(self, material_id: int) -> np.ndarray:
        r""":math:`\Sigma_{2n}` dense ``(ng, ng)`` matrix for one material.

        Cached dense expansion of :attr:`Mixture.Sig2` (sparse upstream).
        """
        if self._n2n_dense is None:
            self._build_dense_caches()
        return self._n2n_dense[material_id]  # type: ignore[index]

    def fission_production_per_material(self, material_id: int) -> np.ndarray:
        r""":math:`\nu \Sigma_f` ``(ng,)`` vector for one material."""
        return self.materials[material_id].SigP

    def chi_per_material(self, material_id: int) -> np.ndarray:
        r""":math:`\chi` ``(ng,)`` vector for one material."""
        return self.materials[material_id].chi

    def _build_dense_caches(self) -> None:
        """Populate dense Legendre scattering + (n,2n) caches.

        Called lazily by :meth:`sig_s_legendre` / :meth:`n2n_matrix`
        on first access.  Caches the dense ``(ng, ng)`` matrices so
        ``apply_p0_in_scatter`` / ``apply_n2n`` / ``apply_legendre_scattering_moments``
        avoid repeated sparse-to-dense conversion in the hot path.
        """
        sig_s_dense: dict[int, list[np.ndarray]] = {}
        n2n_dense: dict[int, np.ndarray] = {}
        for mid, mix in self.materials.items():
            sig_s_dense[mid] = [np.asarray(s.todense()) for s in mix.SigS]
            n2n_dense[mid] = np.asarray(mix.Sig2.todense())
        self._sig_s_dense = sig_s_dense
        self._n2n_dense = n2n_dense

    # ── Convenience metadata ──────────────────────────────────────────

    @property
    def ng(self) -> int:
        """Energy group count — read-through from :attr:`mesh`."""
        return self.mesh.ng

    @property
    def spatial_shape(self) -> tuple[int, ...]:
        """Per-axis cell counts — read-through from :attr:`mesh`.

        C5.2 (#225): replaces the retired ``nx``/``ny`` pair — ONE
        rank-generic metadata read-through instead of two d-bound ones
        (``ny`` lied at d=1 and silently truncated at d≥3).
        """
        return self.mesh.spatial_shape

    # ── Typed verbs — the lifted per-material dispatch loops ─────────

    def apply_p0_in_scatter(self, Q: np.ndarray, phi: np.ndarray) -> None:
        r"""Add P0 in-scatter source :math:`\Sigma_{s,0}^T\,\phi` to ``Q``.

        For each cell ``(ix, iy)`` of material ``mid``,
        :math:`Q[:, ix, iy] \mathrel{+}= \Sigma_{s,0}^{mid,T} \phi[:, ix, iy]`.
        In numpy speak: ``Q[:, ix, iy] += sig_s0[mid].T @ phi[:, ix, iy]``,
        named via :func:`numpy.einsum` to expose the
        source-spectrum-to-sink-spectrum contraction.

        Encapsulates the per-material loop that previously lived at
        ``scattering.py:405`` (``ScatteringOperator.add_iso_source``).

        Spatial-moment-axis-agnostic (#240 D5b-S3 — the
        :math:`\Sigma_s \otimes I_{\rm spatial}` lift): the cell-axis subscript
        carries a trailing ``...`` so a SPATIAL-MOMENT axis (LD's ``2^d`` per
        cell, the diffusion-limit-consistent slope source :math:`\Sigma_s\hat\phi`)
        rides through as a SPECTATOR broadcast — :math:`\Sigma_s` carries NO
        spatial-moment index, so it is applied to EVERY spatial moment of
        :math:`\phi` independently.  At the single-moment closures (DD/Step,
        ``phi`` rank ``(ng, nx, ny)``) the trailing axis is ABSENT and the
        ``...`` matches nothing → BYTE-IDENTICAL to the pre-S3 ``fc->gc``
        (verified rank-2-exact).

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source, shape ``(ng, nx, ny)`` (or ``(ng, nx, ny, 2^d)``
            at an LD multi-moment closure).  Modified in place.
        phi : np.ndarray
            Scalar flux, same shape as ``Q`` (its spatial-moment axis, if any,
            is the spectator).
        """
        for mid, idx in self.cells_by_material.items():
            sig_s0 = self.sig_s_legendre(mid)[0]  # (ng, ng)
            cells = (slice(None), *idx)
            Q[cells] += np.einsum(
                "fg,fc...->gc...", sig_s0, phi[cells],
            )

    def apply_n2n(self, Q: np.ndarray, phi: np.ndarray) -> None:
        r"""Add (n,2n) source :math:`2\,\Sigma_{2n}^T\,\phi` to ``Q``.

        For each cell ``(ix, iy)`` of material ``mid``,
        :math:`Q[:, ix, iy] \mathrel{+}= 2 \cdot \Sigma_{2n}^{mid,T} \phi[:, ix, iy]`.

        Encapsulates the per-material loop that previously lived at
        ``scattering.py:426`` (``ScatteringOperator.add_n2n_source``).

        Spatial-moment-axis-agnostic (#240 D5b-S3, the same
        :math:`\Sigma \otimes I_{\rm spatial}` spectator-broadcast as
        :meth:`apply_p0_in_scatter`): the cell-axis subscript carries a trailing
        ``...`` so an LD spatial-moment axis rides through; byte-identical at the
        single-moment closures (trailing axis absent).

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source, shape ``(ng, nx, ny)`` (or ``(ng, nx, ny, 2^d)``
            at an LD multi-moment closure).  Modified in place.
        phi : np.ndarray
            Scalar flux, same shape as ``Q``.
        """
        for mid, idx in self.cells_by_material.items():
            sig2 = self.n2n_matrix(mid)
            cells = (slice(None), *idx)
            Q[cells] += 2.0 * np.einsum(
                "fg,fc...->gc...", sig2, phi[cells],
            )

    def apply_legendre_scattering_moments(
        self,
        moments: np.ndarray,
        L: int,
        skip_l0: bool,
    ) -> np.ndarray:
        r"""Apply per-ℓ block-diagonal scattering :math:`\Lambda` to a
        moment field.

        Implements

        .. math::

            (\Lambda \phi)_\ell^m(\vec r)\bigg|_g
            = \sum_{g'} \Sigma_{s,\ell}^{m(\vec r)}(g' \to g)\,
                       \phi_\ell^m(\vec r)\bigg|_{g'},

        with the per-material structure folded into the cell axis
        via :attr:`cells_by_material`.  Encapsulates the per-material
        loop that previously lived at ``scattering.py:234``
        (``LegendreMomentScattering.apply``).

        Parameters
        ----------
        moments : np.ndarray
            Moment field of shape ``(L+1, 2L+1, ng, nx, ny)`` (or
            ``(L+1, 2L+1, ng, nx, ny, 2^d)`` at an LD multi-moment closure).  The
            :math:`m` axis is the addition-theorem-shifted index where
            slot ``l + m`` holds :math:`(\ell, m)`.
        L : int
            Maximum Legendre order retained.
        skip_l0 : bool
            When ``True``, leave the :math:`\ell = 0` block as zero —
            the P0 in-scatter goes through :meth:`apply_p0_in_scatter`
            on the reaction-rate fast path.  Set ``False`` when the
            full :math:`R \Lambda M \psi` composition is needed.

        Returns
        -------
        np.ndarray
            Same shape as ``moments``.

        Notes
        -----
        Spatial-moment-axis-agnostic (#240 D5b-S3, the
        :math:`\Sigma_s \otimes I_{\rm spatial}` lift): the cell-axis subscript
        carries a trailing ``...`` so an LD spatial-moment axis rides through as
        a spectator broadcast; byte-identical at the single-moment closures
        (trailing axis absent).
        """
        out = np.zeros_like(moments)
        l_start = 1 if skip_l0 else 0
        for mid, idx in self.cells_by_material.items():
            sig_s_mid = self.sig_s_legendre(mid)
            cells = (Ellipsis, *idx)
            for l in range(l_start, L + 1):
                n_m = 2 * l + 1
                # Trailing-contiguous indexing pattern (see notes in
                # the original LegendreMomentScattering.apply): keeps
                # numpy from rearranging axes when fancy-indexing.
                moments_view = moments[l, :n_m][cells]
                out_block = np.einsum(
                    "mfc...,fg->mgc...", moments_view, sig_s_mid[l],
                )
                out[l, :n_m][cells] = (
                    out_block + out[l, :n_m][cells]
                )
        return out

    def add_n2n_to_group_rate(
        self,
        rate: np.ndarray,
        flux_distribution: np.ndarray,
        volume: np.ndarray,
    ) -> None:
        r"""Add the (n,2n) contribution to a per-group production rate.

        For each material, accumulates
        :math:`2 \int_V \Sigma_{2n}^{m,g'\to g} \phi_{g'}(\vec r)\,dV`
        into ``rate``.  Encapsulates the per-material loop that
        previously lived at ``solver.py:429``
        (``SNSolver.compute_group_production_rate``).

        Parameters
        ----------
        rate : np.ndarray
            Per-group production rate ``(ng,)``.  Modified in place.
        flux_distribution : np.ndarray
            Scalar flux ``(ng, nx, ny)``.
        volume : np.ndarray
            Per-cell volume ``(nx, ny)``.
        """
        for mid, idx in self.cells_by_material.items():
            sig2 = self.n2n_matrix(mid)
            cells = (slice(None), *idx)
            phi_cells_g = flux_distribution[cells].T  # (n_cells, ng)
            n2n_cell_g = 2.0 * (phi_cells_g @ sig2)
            rate += np.einsum("c,cg->g", volume[idx], n2n_cell_g)

    # ── Foldable / residual split (Phase G four-operator algebra) ────
    #
    # The scattering operator's foldable_part / residual_part /
    # is_foldable_into_sigma_r / foldable_sigma methods all required
    # per-material loops over the source-of-truth sig_s dict.  They
    # become typed accessors here, encapsulating the dispatch.

    def foldable_sig_s(self) -> dict[int, list[np.ndarray]]:
        r"""Per-material P0 diagonal-only Legendre lists.

        For each material ``mid``, returns ``[np.diag(np.diag(sig_s[mid][0]))]``
        — the within-group self-scatter cross-section in matrix form.
        Consumed by :meth:`ScatteringOperator.foldable_part` to build
        a sibling operator carrying only diagonal P0.

        Returns
        -------
        dict[int, list[np.ndarray]]
            ``{mid: [diag_only]}`` with one-element lists (P0 only).
        """
        out: dict[int, list[np.ndarray]] = {}
        for mid in self.materials:
            p0 = self.sig_s_legendre(mid)[0]
            out[mid] = [np.diag(np.diag(p0))]
        return out

    def residual_sig_s(self) -> dict[int, list[np.ndarray]]:
        r"""Per-material residual Legendre scattering lists.

        For each material ``mid``, returns
        ``[off_diagonal_P0, Σ_{s,1}, ..., Σ_{s,L}]`` — the cross-group
        P0 (off-diagonal) plus every :math:`\ell \ge 1` block verbatim.
        Consumed by :meth:`ScatteringOperator.residual_part` to build
        a sibling operator carrying the non-foldable channels.

        Returns
        -------
        dict[int, list[np.ndarray]]
            ``{mid: [cross_group_P0, Σ_{s,1}, ...]}``.
        """
        out: dict[int, list[np.ndarray]] = {}
        for mid in self.materials:
            mats = self.sig_s_legendre(mid)
            p0 = mats[0]
            cross_group = p0 - np.diag(np.diag(p0))
            out[mid] = [cross_group, *mats[1:]]
        return out

    def is_p0_diagonal_with_zero_n2n(self) -> bool:
        r"""Return ``True`` iff every material carries diagonal P0 +
        zero (n,2n).

        Structural predicate consumed by
        :meth:`ScatteringOperator.is_foldable_into_sigma_r` (after
        the ``scattering_order == 0`` check).  Uses :func:`numpy.allclose`
        defaults so FP-rounding from ``np.diag(np.diag(...))``
        construction does not produce a false negative.

        Returns
        -------
        bool
            ``True`` iff for every material, ``sig_s[mid][0]`` is
            diagonal AND ``sig2[mid]`` is zero.
        """
        for mid in self.materials:
            p0 = self.sig_s_legendre(mid)[0]
            if not np.allclose(p0, np.diag(np.diag(p0))):
                return False
            if not np.allclose(self.n2n_matrix(mid), 0.0):
                return False
        return True

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""Per-material foldable cross-section :math:`(\sigma_{s,0}^{g\to g})_g`.

        For each material ``mid``, returns the ``(ng,)`` array
        ``np.diag(sig_s[mid][0])`` — the per-group within-group
        self-scatter cross-section that the cell-balance denominator
        absorbs into :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}`.

        Returns
        -------
        dict[int, np.ndarray]
            ``{mid: (ng,) array}``.  Each entry is a fresh copy.
        """
        return {
            mid: np.diag(self.sig_s_legendre(mid)[0]).copy()
            for mid in self.materials
        }

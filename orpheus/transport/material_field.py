r"""Per-material kernel fields — the datum a channel binding retains.

A transport interaction channel has exactly two representation-free
halves: **what each material does** (the channel's kernel — the frozen
per-material datum of :mod:`orpheus.transport.kernels`) and **where each
material sits** (the mesh's material layout,
:attr:`~orpheus.transport.mesh.material_mesh.MaterialMesh.cells_by_material`).
:class:`MaterialField` is that pairing and nothing else: a frozen map
``{material id → kernel}`` plus the shared cell-index partition, with the
ONE gathered per-cell ``(ng, ng)`` contraction primitive every channel
verb rides.

This is CS4c step 3's landing for the posing-filtration **O-6** ruling
(R13): the eight ``MaterialXSField.apply_*`` arms moved here as the
kernels' *array verbs* — the per-material dispatch loop written once, the
einsums ported **verbatim** (bit-identity is gated per arm). Since #426
step 2 (2026-09-04) the scattering and :math:`(n,2n)` channels are two
instances of ONE field, :class:`TransferMaterialField`, and the
multiplicity is read from the kernels' own datum
(:attr:`~orpheus.transport.kernels.TransferKernel.multiplicity`) into
every verb as its ``scale`` — the P0 verbs, the per-ℓ moment verbs and the
group-rate verb alike; a yield of 1 leaves each verb's arithmetic exactly
as it was (the ``scale != 1`` fast path is skipped).

Design notes
------------
* **The layout is handed in, never derived.** ``cells_by_material`` is
  the mesh's own cached partition object, shared by every field over one
  mesh — this class never touches a mesh (no ``np.where``, no
  ``mat_map``), so a field is constructible from pure data in a test
  with a two-line dict.
* **Channel verbs live on channel subclasses.** The base owns the
  pairing, the admission, and the contraction primitive; what an
  interaction *means* (a yield-weighted transfer stack, a rank-1
  emission) is the subclass's vocabulary —
  :class:`TransferMaterialField` and :class:`FissionMaterialField`.
* **Accumulation semantics.** The scalar-carrier verbs ADD into a caller
  accumulator in place (``Q[cells] += …``), exactly like the arms they
  replace — materials partition the cells, so per-material accumulation
  order cannot interact.
* **Spatial-moment spectator.** Every einsum keeps the trailing ``...``
  cell subscript of the arms (#240 D5b-S3): an LD ``2^d`` spatial-moment
  axis rides through as a broadcast spectator; byte-identical when the
  axis is absent.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import TYPE_CHECKING, Generic, Iterator, Mapping, TypeVar

import numpy as np

from orpheus.transport.kernels import FissionKernel, TransferKernel

if TYPE_CHECKING:
    from orpheus.numerics.spaces.moment_head import MomentHead
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = [
    "FissionMaterialField",
    "MaterialField",
    "TransferMaterialField",
]

K = TypeVar("K")

#: The two gathered-contraction spellings every scalar-carrier verb rides.
#: Forward contracts the SOURCE group of a ``[g_from, g_to]`` transfer
#: matrix (output = the sink group ``g``); transpose contracts the SINK
#: group.  The trailing ``...`` is the #240 spatial-moment spectator.
_FORWARD = "fg,fc...->gc..."
_TRANSPOSE = "fg,gc...->fc..."


@dataclass(frozen=True)
class MaterialField(Generic[K]):
    r"""A per-material kernel datum distributed over a mesh's material layout.

    The representation-free pairing a channel binding retains: the frozen
    kernels (one per material id) and the mesh-owned cell partition. Both
    mappings are wrapped read-only at construction; the kernels are frozen
    dataclasses with write-protected arrays, so no consumer can reach
    production data through a field.

    Parameters
    ----------
    per_material : Mapping[int, K]
        The channel kernel of each material, keyed by material id.
    cells_by_material : Mapping[int, tuple[np.ndarray, ...]]
        The mesh's material layout —
        :attr:`MaterialMesh.cells_by_material
        <orpheus.transport.mesh.material_mesh.MaterialMesh.cells_by_material>`
        (one ``np.where`` index tuple per material, arity = mesh ndim).
        Handed in as data; a field never derives it.
    """

    per_material: Mapping[int, K]
    cells_by_material: Mapping[int, tuple[np.ndarray, ...]] = field(repr=False)

    def __post_init__(self) -> None:
        if len(self.per_material) == 0:
            raise ValueError(
                f"{type(self).__name__} needs at least one material; got "
                f"an empty per_material map"
            )
        missing = set(self.cells_by_material) - set(self.per_material)
        if missing:
            raise ValueError(
                f"{type(self).__name__}: the layout names material ids "
                f"{sorted(missing)} that carry no kernel; every laid-out "
                f"material must have one"
            )
        ngs = {k.ng for k in self.per_material.values()}  # type: ignore[attr-defined]
        if len(ngs) != 1:
            raise ValueError(
                f"{type(self).__name__} requires one uniform group "
                f"structure across materials; got ng values {sorted(ngs)}"
            )
        object.__setattr__(
            self, "per_material", MappingProxyType(dict(self.per_material)),
        )
        object.__setattr__(
            self,
            "cells_by_material",
            MappingProxyType(dict(self.cells_by_material)),
        )

    @property
    def ng(self) -> int:
        """Number of energy groups (uniform across materials by admission)."""
        return next(iter(self.per_material.values())).ng  # type: ignore[attr-defined]

    def _laid_out(self) -> Iterator[tuple[K, tuple[np.ndarray, ...]]]:
        """Yield ``(kernel, cell indices)`` per laid-out material — the ONE
        per-material dispatch loop every verb rides."""
        for mid, idx in self.cells_by_material.items():
            yield self.per_material[mid], idx

    def _accumulate_contracted(
        self,
        Q: np.ndarray,
        x: np.ndarray,
        matrix_of,
        *,
        spec: str,
        scale: int = 1,
    ) -> None:
        r"""``Q[cells] += scale · einsum(spec, matrix_of(kernel), x[cells])``
        per laid-out material — the gathered per-cell ``(ng, ng)``
        contraction primitive (the one einsum home of the scalar-carrier
        verbs; the arms' arithmetic, verbatim)."""
        for kernel, idx in self._laid_out():
            cells = (slice(None), *idx)
            contracted = np.einsum(spec, matrix_of(kernel), x[cells])
            if scale != 1:
                contracted = scale * contracted
            Q[cells] += contracted


#: The group contraction of one degree block, by the RANK of the angular head
#: (how many leading axes the moments tensor owes the head): the rectangular
#: harmonic head keeps an ``m`` axis in front of the group axis, the flat
#: Legendre head does not. The letters are the former inline specs, verbatim
#: for the harmonic rows (bit-identical); ``f`` is the source group, ``g`` the
#: sink group, ``c...`` the cells.
_BLOCK_CONTRACTION: dict[tuple[int, bool], str] = {
    (2, False): "mfc...,fg->mgc...",
    (2, True): "mfc...,gf->mgc...",
    (1, False): "fc...,fg->gc...",
    (1, True): "fc...,gf->gc...",
}


def _block_contraction(head: "MomentHead", *, transposed: bool) -> str:
    """The einsum spec for one degree block of a moments tensor with this head."""
    rank = len(head.shape)
    try:
        return _BLOCK_CONTRACTION[(rank, transposed)]
    except KeyError:
        raise NotImplementedError(
            f"moment verbs: an angular head of rank {rank} ({head.name!r}) "
            f"has no block contraction; the shipped heads are rank 2 (the "
            f"real harmonics) and rank 1 (the Legendre basis)."
        ) from None


@dataclass(frozen=True)
class TransferMaterialField(MaterialField[TransferKernel]):
    r"""A transfer channel's field: yield-carrying Legendre stacks over the layout.

    ONE channel per field — every kernel carries the same
    :attr:`multiplicity` (admission), so a field IS the scattering channel
    (:meth:`scattering`, yield 1) or the :math:`(n,2n)` channel
    (:meth:`n2n`, yield :data:`~orpheus.transport.kernels.N2N_MULTIPLICITY`)
    of a facade, and the two differ in that datum alone. The array-verb
    home of the P0 in-transfer and the per-ℓ moment redistribution (the
    former ``MaterialXSField.apply_p0_in_scatter`` / ``apply_n2n`` /
    ``apply_legendre_scattering_moments`` arms and their transposes,
    einsums verbatim) — every verb scales its contraction by the yield.
    The field's own :attr:`order` IS the binding's truncation — a consumer
    that wants :math:`P_L` holds a field :meth:`at_order` :math:`L`, so
    the verb signatures carry no ``L`` parameter (single source; the
    moments tensor must agree, and the moment verbs refuse a mismatch).
    """

    def __post_init__(self) -> None:
        super().__post_init__()
        yields = {k.multiplicity for k in self.per_material.values()}
        if len(yields) != 1:
            raise ValueError(
                f"TransferMaterialField is ONE channel: every material's "
                f"kernel must carry the same multiplicity; got "
                f"{sorted(yields)}"
            )
        # ONE order across materials, by construction: a material whose
        # stack is shorter than the widest is padded to it (exact zeros —
        # the evaluation's own statement; the same law the isotope sum and
        # the binding's at_order apply), so a heterogeneous library (Be-9
        # with NL = 7 beside H-1 with NL = 1) builds a field whose order is
        # readable and whose moment verbs need no lazy refusal.
        widest = max(k.order for k in self.per_material.values())
        object.__setattr__(
            self,
            "per_material",
            MappingProxyType({
                mid: k.at_order(widest) for mid, k in self.per_material.items()
            }),
        )

    @classmethod
    def scattering(cls, mat_xs: "MaterialXSField") -> "TransferMaterialField":
        """The scattering channel of a :class:`MaterialXSField` facade —
        fresh densified kernels (yield 1), the mesh's own layout object."""
        return cls(
            per_material={
                mid: TransferKernel.scattering(mat_xs.materials[mid])
                for mid in mat_xs.materials
            },
            cells_by_material=mat_xs.mesh.cells_by_material,
        )

    @classmethod
    def n2n(cls, mat_xs: "MaterialXSField") -> "TransferMaterialField":
        r"""The :math:`(n,2n)` channel of a :class:`MaterialXSField` facade
        — fresh densified kernels (yield 2), the mesh's own layout object."""
        return cls(
            per_material={
                mid: TransferKernel.n2n(mat_xs.materials[mid])
                for mid in mat_xs.materials
            },
            cells_by_material=mat_xs.mesh.cells_by_material,
        )

    @property
    def order(self) -> int:
        r"""The field's Legendre order :math:`L` — uniform across
        materials by CONSTRUCTION (``__post_init__`` pads every kernel to
        the widest stored order; :meth:`at_order` moves the whole field
        to another)."""
        return next(iter(self.per_material.values())).order

    @property
    def multiplicity(self) -> int:
        r"""The channel's yield :math:`y` (uniform across materials by
        admission) — the ``scale`` of every verb below."""
        return next(iter(self.per_material.values())).multiplicity

    @property
    def is_isotropic(self) -> bool:
        r"""``True`` iff every material's kernel is exactly isotropic
        (:attr:`TransferKernel.is_isotropic
        <orpheus.transport.kernels.TransferKernel.is_isotropic>`) — the
        field's :math:`\Lambda_{\ell\ge 1}` is the zero operator."""
        return all(k.is_isotropic for k in self.per_material.values())

    def at_order(self, order: int) -> "TransferMaterialField":
        r"""The P\ :sub:`order` sub-field — every kernel :meth:`at_order
        <orpheus.transport.kernels.TransferKernel.at_order>` (identity at
        the stored order, truncation below it, exact zeros above it). The
        layout object is shared, not copied."""
        return replace(
            self,
            per_material={
                mid: k.at_order(order)
                for mid, k in self.per_material.items()
            },
        )

    # ── scalar-carrier verbs (the P0 reaction-rate fast path) ─────────

    def add_p0_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        r"""Add the P0 emission :math:`y\,\Sigma_{c,0}^{T}\phi` to ``Q``
        in place (the former ``apply_p0_in_scatter`` / ``apply_n2n`` arms).

        Per cell of material ``m``: ``Q[:, c] += y · p0[m].T @ phi[:, c]``,
        spelled as the source-group contraction of the ``[g_from, g_to]``
        stack head. Shapes ``(ng, *spatial)`` (+ an optional trailing
        spatial-moment spectator axis on both).
        """
        self._accumulate_contracted(
            Q, phi, lambda k: k.p0, spec=_FORWARD, scale=self.multiplicity,
        )

    def add_p0_source_transpose(self, Q: np.ndarray, chi: np.ndarray) -> None:
        r"""Add :math:`y\,\Sigma_{c,0}\,\chi` to ``Q`` in place — the bare
        Euclidean transpose of :meth:`add_p0_source` (sink-group
        contraction; the group-asymmetric factor of the adjoint isotropic
        emission, #276). The metric Hilbert adjoint stays the ``.H``
        wrapper's job."""
        self._accumulate_contracted(
            Q, chi, lambda k: k.p0, spec=_TRANSPOSE, scale=self.multiplicity,
        )

    # ── moment-carrier verbs (the frame-conjugated Λ) ─────────────────

    def _moment_blocks(
        self, moments: np.ndarray, *, skip_l0: bool, head: "MomentHead", spec: str,
    ) -> np.ndarray:
        r"""The per-ℓ block-diagonal group contraction shared by the two
        moment verbs (the former ``apply_legendre_scattering_moments``
        loop — trailing-contiguous indexing kept so numpy never
        rearranges axes under the fancy index), every block scaled by the
        channel's yield.

        The moments tensor's leading axes are the angular HEAD's
        (:class:`~orpheus.numerics.spaces.moment_head.MomentHead`): the
        rectangular ``(L+1, 2L+1)`` of the real harmonics, or the flat
        ``(L+1,)`` of the Legendre basis on a 1-D rule (#429, 2026-09-02).
        The head says which index tuple is the degree-:math:`\ell` block
        and how many axes precede the group axis; until 2026-09-02 this
        loop spelled the harmonics' m-axis into its einsum and its
        slicing, and a flat head would have contracted the GROUP axis as
        if it were :math:`m` — silently.
        """
        L = self.order
        rank = len(head.shape)
        ng = self.ng
        if (
            moments.ndim < rank + 1
            or moments.shape[:rank] != tuple(head.shape)
            or moments.shape[rank] != ng
        ):
            raise ValueError(
                f"TransferMaterialField moment verb: the moments tensor "
                f"{moments.shape} must lead with the angular head "
                f"{tuple(head.shape)} of {head.name!r} followed by the "
                f"{ng}-group axis (this field is at order {L}); bring the "
                f"field to the operator's order and hand the operator's own "
                f"head — a rectangular tensor offered with a flat head would "
                f"otherwise read its m-axis as the groups."
            )
        scale = self.multiplicity
        out = np.zeros_like(moments)
        l_start = 1 if skip_l0 else 0
        for kernel, idx in self._laid_out():
            cells = (slice(None),) * rank + tuple(idx)
            for l in range(l_start, L + 1):
                block = head.degree_block(l)
                moments_view = moments[block][cells]
                out_block = np.einsum(spec, moments_view, kernel.moments[l])
                if scale != 1:
                    out_block = scale * out_block
                out[block][cells] = out_block + out[block][cells]
        return out

    def moment_source(
        self, moments: np.ndarray, *, skip_l0: bool, head: "MomentHead",
    ) -> np.ndarray:
        r"""Apply the per-ℓ block-diagonal redistribution :math:`\Lambda`
        to a moment field (the former ``apply_legendre_scattering_moments``
        / ``apply_n2n_moments`` arms):
        :math:`(\Lambda\phi)_\ell^m|_g = y \sum_{g'}
        \Sigma_{c,\ell}^{m(\vec r)}(g'\to g)\,\phi_\ell^m|_{g'}`.
        ``skip_l0=True`` leaves the :math:`\ell=0` block zero (the P0
        fast path owns it). ``head`` is the operator's angular head (its
        domain), which fixes the layout the contraction reads."""
        return self._moment_blocks(
            moments, skip_l0=skip_l0, head=head,
            spec=_block_contraction(head, transposed=False),
        )

    def moment_source_transpose(
        self, moments: np.ndarray, *, skip_l0: bool, head: "MomentHead",
    ) -> np.ndarray:
        r"""Apply :math:`\Lambda^{T}` — the per-ℓ group-transpose twin of
        :meth:`moment_source` (sink-group contraction ``gf``)."""
        return self._moment_blocks(
            moments, skip_l0=skip_l0, head=head,
            spec=_block_contraction(head, transposed=True),
        )

    # ── the group-rate verb (the k-balance accounting) ────────────────

    def add_to_group_rate(
        self,
        rate: np.ndarray,
        flux_distribution: np.ndarray,
        volume: np.ndarray,
    ) -> None:
        r"""Add the channel's volume-integrated P0 emission to a per-group
        rate (the former ``add_n2n_to_group_rate`` arm): per material,
        :math:`y\int_V \Sigma_{c,0}^{g'\to g}\,\phi_{g'}\,dV` accumulated
        into ``rate`` (shape ``(ng,)``, in place). For :math:`(n,2n)` this
        is the emission the k-balance's net removal subtracts and the
        ERR-052 scale anchor adds; for scattering it is the group
        in-scatter rate."""
        scale = self.multiplicity
        for kernel, idx in self._laid_out():
            cells = (slice(None), *idx)
            phi_cells_g = flux_distribution[cells].T  # (n_cells, ng)
            cell_g = phi_cells_g @ kernel.p0
            if scale != 1:
                cell_g = scale * cell_g
            rate += np.einsum("c,cg->g", volume[idx], cell_g)


@dataclass(frozen=True)
class FissionMaterialField(MaterialField[FissionKernel]):
    r"""The fission channel's field: validated factor pairs over the layout.

    :class:`~orpheus.transport.kernels.FissionKernel`'s first production
    consumer (CS4c step 4): each material's :math:`(\chi, \nu\Sigma_f)`
    pair enters HERE through the kernel's own constructor, so the
    χ simplex/null law runs per material by construction (Pattern 4 —
    a field holding an invalid spectrum is not a value that exists).

    Unlike its transfer sibling this field carries **gather verbs, not
    accumulation verbs**: fission's arithmetic home is the rank-1 dyad
    :math:`|\chi\rangle\langle\nu\Sigma_f|` realized once on the CELLWISE
    factors (the :class:`~orpheus.numerics.operator.RankOneOperator`
    route the bound operators cache at construction), so what the binding
    needs from the datum is the densified factor — bit-identical to the
    facade's per-cell views, being the same pure index gather of the same
    per-material vectors. A per-material contraction verb here would be
    a SECOND spelling of the dyad arithmetic (Pattern 2).
    """

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField",
    ) -> "FissionMaterialField":
        """Extract the fission channel of a :class:`MaterialXSField`
        facade — fresh validated kernels, the mesh's own layout object."""
        return cls(
            per_material={
                mid: FissionKernel.from_mixture(mat_xs.materials[mid])
                for mid in mat_xs.materials
            },
            cells_by_material=mat_xs.mesh.cells_by_material,
        )

    def gather_chi(self, spatial_shape: tuple[int, ...]) -> np.ndarray:
        r"""The cellwise emission spectrum :math:`\chi(\vec r, g)` —
        shape ``(ng, *spatial_shape)``, write-protected.

        Each material's ``(ng,)`` spectrum broadcast into its cells (a
        pure index gather — no arithmetic, so the values are
        bit-identical to any other densification of the same
        per-material vectors)."""
        return self._gather_vector(lambda k: k.chi, spatial_shape)

    def gather_nu_sig_f(self, spatial_shape: tuple[int, ...]) -> np.ndarray:
        r"""The cellwise production cross section
        :math:`\nu\Sigma_f(\vec r, g)` — shape ``(ng, *spatial_shape)``,
        write-protected (see :meth:`gather_chi`)."""
        return self._gather_vector(lambda k: k.nu_sig_f, spatial_shape)

    def _gather_vector(
        self, vector_of, spatial_shape: tuple[int, ...],
    ) -> np.ndarray:
        """Densify a per-material ``(ng,)`` vector over the layout.

        The shape is HANDED IN (by the binding, from its own space's
        bulk shape) — the base-class rule that a field never touches a
        mesh. Refuses a layout index that falls outside the shape, so a
        wrong-space gather fails loudly instead of silently truncating.
        """
        out = np.zeros((self.ng, *spatial_shape), dtype=float)
        for kernel, idx in self._laid_out():
            if len(idx) != len(spatial_shape):
                raise ValueError(
                    f"FissionMaterialField gather: the layout indexes "
                    f"{len(idx)} spatial axes but the requested shape "
                    f"{spatial_shape} has {len(spatial_shape)}"
                )
            out[(slice(None), *idx)] = vector_of(kernel)[:, None]
        out.setflags(write=False)
        return out

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
(R13): the eight ``MaterialXSField.apply_*`` arms move here as the
kernels' *array verbs* — the per-material dispatch loop written once, the
einsums ported **verbatim** (bit-identity is gated per arm), and the
:math:`(n,2n)` multiplicity read from its one home
(:attr:`~orpheus.transport.kernels.N2NKernel.multiplicity`) instead of an
inline literal (XD-2).

Design notes
------------
* **The layout is handed in, never derived.** ``cells_by_material`` is
  the mesh's own cached partition object, shared by every field over one
  mesh — this class never touches a mesh (no ``np.where``, no
  ``mat_map``), so a field is constructible from pure data in a test
  with a two-line dict.
* **Channel verbs live on channel subclasses.** The base owns the
  pairing, the admission, and the contraction primitive; what an
  interaction *means* (P0 in-scatter, Legendre moment redistribution,
  multiplicity-weighted emission) is the subclass's vocabulary —
  :class:`ScatteringMaterialField`, :class:`N2NMaterialField`, and
  :class:`FissionMaterialField` (step 4).
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

from orpheus.transport.kernels import (
    FissionKernel,
    N2NKernel,
    ScatteringKernel,
)

if TYPE_CHECKING:
    from orpheus.numerics.spaces.moment_head import MomentHead
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = [
    "FissionMaterialField",
    "MaterialField",
    "N2NMaterialField",
    "ScatteringMaterialField",
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
        scale: float = 1.0,
    ) -> None:
        r"""``Q[cells] += scale · einsum(spec, matrix_of(kernel), x[cells])``
        per laid-out material — the gathered per-cell ``(ng, ng)``
        contraction primitive (the one einsum home of the scalar-carrier
        verbs; the arms' arithmetic, verbatim)."""
        for kernel, idx in self._laid_out():
            cells = (slice(None), *idx)
            contracted = np.einsum(spec, matrix_of(kernel), x[cells])
            if scale != 1.0:
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
class ScatteringMaterialField(MaterialField[ScatteringKernel]):
    r"""The scattering channel's field: Legendre transfer stacks over the layout.

    The array-verb home of the P0 in-scatter and the per-ℓ moment
    redistribution (the former ``MaterialXSField.apply_p0_in_scatter`` /
    ``…_transpose`` / ``apply_legendre_scattering_moments`` /
    ``…_transpose`` arms, einsums verbatim). The field's own
    :attr:`order` IS the binding's truncation — a consumer that wants
    :math:`P_L` holds a field :meth:`truncated` to :math:`L`, so the verb
    signatures carry no ``L`` parameter (single source; the moments
    tensor must agree, and the moment verbs refuse a mismatch).
    """

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField",
    ) -> "ScatteringMaterialField":
        """Extract the scattering channel of a :class:`MaterialXSField`
        facade — fresh densified kernels, the mesh's own layout object."""
        return cls(
            per_material={
                mid: ScatteringKernel.from_mixture(mat_xs.materials[mid])
                for mid in mat_xs.materials
            },
            cells_by_material=mat_xs.mesh.cells_by_material,
        )

    @property
    def order(self) -> int:
        r"""The field's Legendre truncation :math:`L` (uniform across
        materials by admission)."""
        orders = {k.order for k in self.per_material.values()}
        if len(orders) != 1:
            raise ValueError(
                f"ScatteringMaterialField carries mixed truncation orders "
                f"{sorted(orders)}; truncate to a uniform order first"
            )
        return orders.pop()

    def truncated(self, order: int) -> "ScatteringMaterialField":
        r"""The P\ :sub:`order` sub-field — every kernel truncated
        (:meth:`ScatteringKernel.truncated
        <orpheus.transport.kernels.ScatteringKernel.truncated>` semantics:
        identity at the stored order, refusal above it). The layout object
        is shared, not copied."""
        return replace(
            self,
            per_material={
                mid: k.truncated(order)
                for mid, k in self.per_material.items()
            },
        )

    # ── scalar-carrier verbs (the P0 reaction-rate fast path) ─────────

    def add_p0_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        r"""Add the P0 in-scatter source :math:`\Sigma_{s,0}^{T}\phi` to
        ``Q`` in place (the former ``apply_p0_in_scatter`` arm).

        Per cell of material ``m``: ``Q[:, c] += p0[m].T @ phi[:, c]``,
        spelled as the source-group contraction of the ``[g_from, g_to]``
        stack head. Shapes ``(ng, *spatial)`` (+ an optional trailing
        spatial-moment spectator axis on both).
        """
        self._accumulate_contracted(
            Q, phi, lambda k: k.p0, spec=_FORWARD,
        )

    def add_p0_source_transpose(self, Q: np.ndarray, chi: np.ndarray) -> None:
        r"""Add :math:`\Sigma_{s,0}\,\chi` to ``Q`` in place — the bare
        Euclidean transpose of :meth:`add_p0_source` (sink-group
        contraction; the group-asymmetric factor of the adjoint isotropic
        in-scatter, #276). The metric Hilbert adjoint stays the
        ``.H`` wrapper's job."""
        self._accumulate_contracted(
            Q, chi, lambda k: k.p0, spec=_TRANSPOSE,
        )

    # ── moment-carrier verbs (the frame-conjugated Λ) ─────────────────

    def _moment_blocks(
        self, moments: np.ndarray, *, skip_l0: bool, head: "MomentHead", spec: str,
    ) -> np.ndarray:
        r"""The per-ℓ block-diagonal group contraction shared by the two
        moment verbs (the former ``apply_legendre_scattering_moments``
        loop — trailing-contiguous indexing kept so numpy never
        rearranges axes under the fancy index).

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
        ng = next(iter(self.per_material.values())).ng
        if (
            moments.ndim < rank + 1
            or moments.shape[:rank] != tuple(head.shape)
            or moments.shape[rank] != ng
        ):
            raise ValueError(
                f"ScatteringMaterialField moment verb: the moments tensor "
                f"{moments.shape} must lead with the angular head "
                f"{tuple(head.shape)} of {head.name!r} followed by the "
                f"{ng}-group axis (this field is truncated to order {L}); "
                f"truncate the field to the operator's order and hand the "
                f"operator's own head — a rectangular tensor offered with a "
                f"flat head would otherwise read its m-axis as the groups."
            )
        out = np.zeros_like(moments)
        l_start = 1 if skip_l0 else 0
        for kernel, idx in self._laid_out():
            cells = (slice(None),) * rank + tuple(idx)
            for l in range(l_start, L + 1):
                block = head.degree_block(l)
                moments_view = moments[block][cells]
                out_block = np.einsum(spec, moments_view, kernel.moments[l])
                out[block][cells] = out_block + out[block][cells]
        return out

    def moment_source(
        self, moments: np.ndarray, *, skip_l0: bool, head: "MomentHead",
    ) -> np.ndarray:
        r"""Apply the per-ℓ block-diagonal redistribution :math:`\Lambda`
        to a moment field (the former
        ``apply_legendre_scattering_moments`` arm):
        :math:`(\Lambda\phi)_\ell^m|_g = \sum_{g'}
        \Sigma_{s,\ell}^{m(\vec r)}(g'\to g)\,\phi_\ell^m|_{g'}`.
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


@dataclass(frozen=True)
class N2NMaterialField(MaterialField[N2NKernel]):
    r"""The :math:`(n,2n)` channel's field: reaction matrices over the layout.

    The array-verb home of the multiplicity-weighted emission (the former
    ``MaterialXSField.apply_n2n`` / ``…_transpose`` /
    ``apply_n2n_moments`` / ``…_transpose`` /
    ``add_n2n_to_group_rate`` arms, einsums verbatim). The multiplicity
    enters every emission verb from its ONE home —
    :attr:`N2NKernel.multiplicity
    <orpheus.transport.kernels.N2NKernel.multiplicity>` — never as an
    inline literal (XD-2): the kernel stores the raw reaction matrix
    :math:`\Sigma_{2n}` and the verbs weight it, exactly as the
    loss/emission channel ruling requires.
    """

    @classmethod
    def from_material_xs(cls, mat_xs: "MaterialXSField") -> "N2NMaterialField":
        """Extract the :math:`(n,2n)` channel of a :class:`MaterialXSField`
        facade — fresh densified kernels, the mesh's own layout object."""
        return cls(
            per_material={
                mid: N2NKernel.from_mixture(mat_xs.materials[mid])
                for mid in mat_xs.materials
            },
            cells_by_material=mat_xs.mesh.cells_by_material,
        )

    def add_emission(self, Q: np.ndarray, phi: np.ndarray) -> None:
        r"""Add the :math:`(n,2n)` source
        :math:`\nu_{2n}\,\Sigma_{2n}^{T}\phi` to ``Q`` in place (the
        former ``apply_n2n`` arm; :math:`\nu_{2n}` =
        :attr:`N2NKernel.multiplicity
        <orpheus.transport.kernels.N2NKernel.multiplicity>`)."""
        self._accumulate_contracted(
            Q, phi, lambda k: k.matrix,
            spec=_FORWARD, scale=float(N2NKernel.multiplicity),
        )

    def add_emission_transpose(self, Q: np.ndarray, chi: np.ndarray) -> None:
        r"""Add :math:`\nu_{2n}\,\Sigma_{2n}\,\chi` to ``Q`` in place —
        the bare Euclidean transpose of :meth:`add_emission`."""
        self._accumulate_contracted(
            Q, chi, lambda k: k.matrix,
            spec=_TRANSPOSE, scale=float(N2NKernel.multiplicity),
        )

    def moment_emission(self, moments: np.ndarray, *, head: "MomentHead") -> np.ndarray:
        r"""Apply the :math:`\ell=0` moment operator
        :math:`\nu_{2n}\,\Sigma_{2n}` (the former ``apply_n2n_moments``
        arm): only the ``[0, 0]`` block is read and written; every
        :math:`\ell\ge 1` block stays zero — because this field's
        kernels carry ONE matrix
        (:meth:`~orpheus.transport.kernels.N2NKernel.from_mixture`
        selects ``Sig2[0]``), ORPHEUS modelling the :math:`(n,2n)`
        emission at :math:`P_0`.  ⚠ That is a modelling choice, NOT a
        property of the reaction — and since #426 step 1 it is an
        OPERATOR-layer choice, not a data-layer loss: ``Mixture.Sig2``
        carries every Legendre order the tape stores, so the moments
        this verb leaves at zero are available one call upstream
        (``docs/theory/methods/sn/adjoint.rst`` §sn-n2n-p0-truncation,
        issue #426 step 2)."""
        return self._moment_l0(
            moments, head=head, spec=_block_contraction(head, transposed=False),
        )

    def moment_emission_transpose(
        self, moments: np.ndarray, *, head: "MomentHead",
    ) -> np.ndarray:
        r"""Apply :math:`(\nu_{2n}\,\Sigma_{2n})^{T}` — the ℓ=0
        group-transpose twin of :meth:`moment_emission`."""
        return self._moment_l0(
            moments, head=head, spec=_block_contraction(head, transposed=True),
        )

    def _moment_l0(
        self, moments: np.ndarray, *, head: "MomentHead", spec: str,
    ) -> np.ndarray:
        # Only the l = 0 block is read and written; its index tuple and the
        # number of axes before the group axis are the HEAD's to say
        # (rectangular harmonics or the flat Legendre head — #429).
        mult = float(N2NKernel.multiplicity)
        rank = len(head.shape)
        block = head.degree_block(0)
        out = np.zeros_like(moments)
        for kernel, idx in self._laid_out():
            cells = (slice(None),) * rank + tuple(idx)
            mv = moments[block][cells]  # the l = 0 moment over these cells
            out[block][cells] = (
                mult * np.einsum(spec, mv, kernel.matrix) + out[block][cells]
            )
        return out

    def add_to_group_rate(
        self,
        rate: np.ndarray,
        flux_distribution: np.ndarray,
        volume: np.ndarray,
    ) -> None:
        r"""Add the :math:`(n,2n)` contribution to a per-group production
        rate (the former ``add_n2n_to_group_rate`` arm): per material,
        :math:`\nu_{2n}\int_V \Sigma_{2n}^{g'\to g}\,\phi_{g'}\,dV`
        accumulated into ``rate`` (shape ``(ng,)``, in place)."""
        mult = float(N2NKernel.multiplicity)
        for kernel, idx in self._laid_out():
            cells = (slice(None), *idx)
            phi_cells_g = flux_distribution[cells].T  # (n_cells, ng)
            n2n_cell_g = mult * (phi_cells_g @ kernel.matrix)
            rate += np.einsum("c,cg->g", volume[idx], n2n_cell_g)


@dataclass(frozen=True)
class FissionMaterialField(MaterialField[FissionKernel]):
    r"""The fission channel's field: validated factor pairs over the layout.

    :class:`~orpheus.transport.kernels.FissionKernel`'s first production
    consumer (CS4c step 4): each material's :math:`(\chi, \nu\Sigma_f)`
    pair enters HERE through the kernel's own constructor, so the
    χ simplex/null law runs per material by construction (Pattern 4 —
    a field holding an invalid spectrum is not a value that exists).

    Unlike its scattering/(n,2n) siblings this field carries **gather
    verbs, not accumulation verbs**: fission's arithmetic home is the
    rank-1 dyad :math:`|\chi\rangle\langle\nu\Sigma_f|` realized once on
    the CELLWISE factors (the
    :class:`~orpheus.numerics.operator.RankOneOperator` route the bound
    operators cache at construction), so what the binding needs from the
    datum is the densified factor — bit-identical to the facade's
    per-cell views, being the same pure index gather of the same
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

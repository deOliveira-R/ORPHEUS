r"""The model-independent isotropic energy operators on the scalar flux.

Every transport model — discrete ordinates (SN), collision probability (CP),
method of characteristics (MoC), diffusion, the homogeneous 0-D solver, Monte
Carlo — builds the **isotropic** (:math:`\ell=0`) emission source by the SAME
per-cell energy operation on the **scalar flux** :math:`\phi`:

.. math::

    Q^{\rm iso}_g(\vec r) \;=\; \sum_{g'} \Sigma_{s,0}(g'\to g)\,\phi_{g'}(\vec r)
        \;+\; 2\,\sum_{g'} \Sigma_{2n}(g'\to g)\,\phi_{g'}(\vec r) .

What differs per model is ONLY the *angular wrapper* around it — how :math:`\phi`
is obtained from the model's flux (SN: :math:`\phi=\sum_n w_n\psi_n`; diffusion /
homogeneous: :math:`\phi` is the native unknown; CP / MoC: per region / per FSR)
and how the scalar source is embedded back (SN: :math:`/W` isotropic broadcast;
diffusion / homogeneous: fold :math:`-\Sigma_{s,0}^{T}` into the loss matrix
:math:`A`; CP / MoC: through the transport kernel). The **energy operation is
model-independent** (cross-domain-attacker + explorer, 2026-06-28: the same
``einsum("fg,f->g")`` was re-implemented 6× across the solver families — a
Cardinal-Rule-2 single-concept-many-places situation).

This module owns that shared energy operation as ONE
:class:`~orpheus.numerics.operator.LinearOperator` on the scalar flux —
:class:`IsotropicTransfer`, the P0 energy binding of a transfer channel
:math:`y\,\Sigma_{c,0}^{T}` — with the two terms of the algebra as its role
subclasses (#426 step 2, 2026-09-04; the F3 ruling: the kernel tier names the
mathematical object, the operator tier names the TERM):

* :class:`IsotropicScattering` — :math:`\Sigma_{s,0}` (P0 in-scatter, yield 1);
* :class:`IsotropicN2N` — :math:`2\,\Sigma_{2n}` (the (n,2n) doubling — a
  *multiplication* channel whose yield also feeds the keff production
  accounting, so it stays its own TERM, summed with
  :class:`IsotropicScattering` for the isotropic emission source).

The roles carry NOTHING but the extraction classmethod (which ``Mixture``
channel the tier-2 mint reads) and the role name — an AST gate
(``tests/transport/test_transfer_roles.py``) asserts it, so the
member-for-member twins these two classes were until 2026-09-04 cannot regrow.

Both are the **scalar (:math:`\ell=0`) realization** of the moment-space operator
:class:`~orpheus.transport.operators.transfer.LegendreMomentTransfer` (at
:math:`\ell=0`): they route through the SAME per-material kernel-field verb
(:meth:`~orpheus.transport.material_field.TransferMaterialField.add_p0_source`
+ its ``…_transpose`` twin — CS4c step 3, the O-6 landing), so the cross-section
DATA, the per-material dispatch and the yield live ONCE (``coding-elegance``
Pattern 2). The harmonic-moment
:attr:`~orpheus.transport.operators.transfer.TransferOperator.full_transfer_kernel`
is the permanent verification oracle for this scalar form.

Layout & carriers
=================

Bare ``np.ndarray`` in / out — the model-portable contract (CP / MoC / diffusion
feed raw scalar-flux arrays; SN passes ``phi.values`` and re-wraps). The action is
**spatial-moment-axis-agnostic** (the trailing ``…`` rides an LD ``2^d`` spectator
axis through, #240 D5b-S3), so a φ̂-carrying scalar flux scatters every spatial
moment by the same :math:`\Sigma_c`.

Capabilities
============

apply + a working ``apply_transpose``. No inverse — the per-cell group-transfer
matrix is generally singular (a thermal group with no up-scatter source has a zero
row); it is *applied*, never inverted. The transpose IS advertised (campaign #276):
:math:`y\,\Sigma_{c,0}^{T}` (the group-axis flip) is the group-asymmetric factor of
the adjoint isotropic emission source.

References
----------

* explorer (2026-06-28): cross-model grounding — the 6× / 5× duplication inventory.
* cross-domain-attacker (2026-06-28):
  ``iso_source_frame_conjugation_unification.md`` — relocate the shared energy
  operator; do NOT mint a ``ConstantBasis`` iso-frame (it forks ``R∘M``).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    LinearOperator,
)
from orpheus.transport.material_field import (
    FissionMaterialField,
    TransferMaterialField,
)
from orpheus.transport.operators.bound_operator import BoundOperator
# Runtime import for the composite-arm isinstance parse (mirrors
# fission.py / multiplication_operator.py): ``FullField`` is a leaf in
# the transport dependency graph (it imports no operators), so this
# module-level import is cycle-free.
from orpheus.transport.full_field import FullField

if TYPE_CHECKING:
    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.operator import TensorProductOperator
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.transport.reaction_rate_functional import ReactionRateFunctional


__all__ = [
    "IsotropicTransfer",
    "IsotropicScattering",
    "IsotropicN2N",
    "IsotropicFission",
]


def _values_of(phi: "np.ndarray | object") -> np.ndarray:
    """Read the bare ``(ng, *spatial[, 2^d])`` array off a flux carrier or ndarray."""
    return np.asarray(getattr(phi, "values", phi))


def _scalar_composite_source(op: "LinearOperator", psi: FullField) -> FullField:
    r"""The shared scalar-composite arm of the iso energy operators (#290 P4).

    Parses the composite's bulk as the SCALAR family (the iso operators
    are scalar-flux operators — an angular composite routes through
    :class:`~orpheus.transport.operators.transfer.TransferOperator`,
    never here), runs the operator's own bare-ndarray kernel on the bulk
    values (single source of truth — the SAME per-cell energy matmul),
    and wraps the result as the closed scalar source composite:
    :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
    bulk ⊕ implicit-zero
    :class:`~orpheus.transport.source_sinks.scalar_boundary_source_sink.ScalarBoundarySourceSink`
    boundary (an energy-transfer operator has no face-trace action —
    the emission term is a CELL quantity). Shared by
    :class:`IsotropicTransfer` and :class:`IsotropicFission` (Pattern 2:
    one arm body, two kernels through ``op.apply``).
    """
    from orpheus.transport.fields._bases import ScalarField
    from orpheus.transport.source_sinks import (
        ScalarBoundarySourceSink,
        ScalarSourceSink,
    )

    bulk = psi.interior
    if not isinstance(bulk, ScalarField):
        raise TypeError(
            f"{type(op).__name__} composite apply: scalar-family bulk "
            f"required (the iso energy operators act on the scalar flux; "
            f"an angular composite's emission routes through "
            f"TransferOperator); got {type(bulk).__name__}."
        )
    # CS4b S4 — the space route: output blocks ride the operand's blocks.
    return FullField(
        interior=ScalarSourceSink(values=op.apply(bulk.values), space=bulk.space),
        boundary=ScalarBoundarySourceSink.zeros(psi.boundary.space),
    )


def _iso_is_assemblable(op: "IsotropicTransfer") -> bool:
    r"""The assembly-axis predicate of the iso energy operator.

    Emittable iff the composite flat layout is known AND the declared
    bulk is the SCALAR family's ``(ng, *cells)`` energy-first tensor —
    the same contract the bare-ndarray kernel arm assumes (a
    block-bearing :class:`FullFieldSpace` whose bulk leading axis is
    the group axis). A scalar-space binding (the diffusion / SN
    internal consumers — ends mandatory since CS4c step 3) honestly
    refuses: no composite flat layout to emit into.
    """
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    space = op.domain
    return (
        isinstance(space, FullFieldSpace)
        and space.interior_space is not None
        and len(space.interior_space.shape) >= 1
        and int(space.interior_space.shape[0]) == int(op.data_ng)
    )


def _assemble_iso_energy_operator(
    op: "IsotropicTransfer",
) -> "SparseAssembledOperator":
    r"""The assembly arm of the iso energy operator (the
    :func:`_scalar_composite_source` sibling).

    **Group-impulse extraction — one source with the production
    kernel.** The per-cell energy blocks are read THROUGH the
    operator's own bare-ndarray ``apply`` (the einsum kernel
    ``add_p0_source``): impulse ``g'`` is the field ``e_{g'} ⊗
    1_cells``, whose image column IS every cell's block column
    ``M_cell[:, g']`` exactly (unit inputs make the kernel's products
    coefficient reads — no summation error). ng kernel calls emit the
    whole cell-block-diagonal bulk. Deliberately NOT
    ``dense_per_material()`` — that is the storage-side
    transpose-convention ORACLE (vv L11) and must stay
    realization-independent of production consumption.

    The emission is bulk-only (an energy-transfer operator has no
    face-trace action — the same zero-boundary fact the composite arm
    encodes), scattered into the composite flat layout ``[bulk C-ravel
    | trace]``.
    """
    from scipy import sparse

    from orpheus.numerics.assembled_operator import SparseAssembledOperator
    from orpheus.numerics.operator import MissingAssembly
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    space = op.domain
    if not _iso_is_assemblable(op):
        raise MissingAssembly(
            f"{type(op).__name__}.assemble requires a block-bearing "
            f"FullFieldSpace with the scalar (ng, *cells) energy-first "
            f"bulk; this instance is space-anonymous or its bulk leading "
            f"axis is not the group axis."
        )
    assert isinstance(space, FullFieldSpace)  # narrowed by the predicate
    assert space.interior_space is not None
    interior_shape = tuple(space.interior_space.shape)
    ng = interior_shape[0]
    n_cells = int(np.prod(interior_shape[1:])) if len(interior_shape) > 1 else 1
    n_interior = ng * n_cells
    n_total = int(space.shape[0])

    bulk_rows = np.arange(n_interior)
    cell_of_row = bulk_rows % n_cells
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    for g_from in range(ng):
        impulse = np.zeros(interior_shape)
        impulse[g_from] = 1.0
        image = np.asarray(op.apply(impulse)).ravel()   # (n_interior,) g-major
        nonzero = image != 0.0
        rows.append(bulk_rows[nonzero])
        cols.append((g_from * n_cells + cell_of_row)[nonzero])
        vals.append(image[nonzero])

    matrix = sparse.coo_array(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_total, n_total),
    )
    return SparseAssembledOperator(matrix, domain=space, codomain=space)


@dataclass(eq=False)
class IsotropicTransfer(BoundOperator):
    r"""The P0 energy binding of a transfer channel: :math:`y\,\Sigma_{c,0}^{T}` on the scalar flux.

    Per cell of material :math:`m`, :math:`(y\,\Sigma_{c,0}^{T}\phi)_g =
    y\sum_{g'}\Sigma_{c,0}^m(g'\to g)\,\phi_{g'}` — the per-material
    group-transfer matmul, vectorised over cells, scaled by the channel's
    yield. The model-independent isotropic emission of ONE channel; the
    two terms of the algebra are the role subclasses
    :class:`IsotropicScattering` (:math:`y = 1`) and :class:`IsotropicN2N`
    (:math:`y = 2`), which add nothing but which channel they read. The
    per-material dispatch and the yield live on the field's verbs
    (:meth:`~orpheus.transport.material_field.TransferMaterialField.add_p0_source`).

    Parameters
    ----------
    transfer : TransferMaterialField
        The channel's kernel field — the per-material Legendre stacks with
        their yield over the mesh layout (only the ``p0`` head is consumed;
        the tier-2 mints bring the field to order 0 so the instance
        retains exactly what it reads).
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once —
        :class:`~orpheus.transport.operators.bound_operator.BoundOperator`,
        CS4c step 3): the scalar-flux space, both — the endomorphism
        sugar lives on the roles' ``from_material_xs``. The OperatorSum
        composition guard and the per-END ng-conformity admission both
        read them; the pre-CS4c ``space=None`` anonymity is retired.
    """

    transfer: "TransferMaterialField"
    # A BULK energy operator (the scalar flux is the bulk block); no boundary action.
    # Class-level constant (unannotated so the dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # CS4a K2 → CS4c per-END form: refuse ends whose EnergyAxis
        # contradicts the data's ng (reach + declared inertness:
        # _energy_conformity docstring).
        self._assert_energy_extent_both_ends(
            self.transfer.ng, operator=type(self).__name__,
        )

    @property
    def data_ng(self) -> int:
        """The bound data's group count (the assembly helpers' read)."""
        return self.transfer.ng

    @property
    def is_adjointable(self) -> bool:
        # apply_transpose realises y·Σ_{c,0}χ (the group-flip A^T); caps ⊇
        # apply_transpose. is_invertible inherits base False.
        return True

    def apply(self, phi: "np.ndarray | FullField | object") -> "np.ndarray | FullField":
        r""":math:`y\,\Sigma_{c,0}^{T}\phi` — the per-cell P0 emission.

        Bare ndarray / flux-carrier in → bare ndarray out (the
        model-portable contract). A scalar :class:`FullField` composite
        in → the closed scalar source composite out (#290 P4 — see
        :func:`_scalar_composite_source`; the kernel is the SAME bare
        arm either way).
        """
        if isinstance(phi, FullField):
            return _scalar_composite_source(self, phi)
        arr = _values_of(phi)
        out = np.zeros_like(arr)
        self.transfer.add_p0_source(out, arr)
        return out

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r""":math:`y\,\Sigma_{c,0}\chi` — the group-flip transpose (the bare Euclidean :math:`A^{T}`).

        Bare-ndarray surface only: the composite transpose arm lands
        with its first consumer (the adjoint diffusion chain, #281) —
        refuse a composite loudly rather than mis-reading it.
        """
        if isinstance(chi, FullField):
            raise TypeError(
                f"{type(self).__name__}.apply_transpose: composite "
                f"FullField transpose is not yet wired (lands with the "
                f"adjoint diffusion consumer, #281); pass the bulk "
                f"values."
            )
        arr = _values_of(chi)
        out = np.zeros_like(arr)
        self.transfer.add_p0_source_transpose(out, arr)
        return out

    # ── The assembly mode (stencil-assembly 2b) ────────────────────────

    @property
    def is_assemblable(self) -> bool:
        # Shared predicate — see :func:`_iso_is_assemblable`.
        return _iso_is_assemblable(self)

    def assemble(self) -> "SparseAssembledOperator":
        r""":math:`[y\,\Sigma_{c,0}^{T}]` — cell-block-diagonal bulk emission
        through the production kernel (:func:`_assemble_iso_energy_operator`)."""
        return _assemble_iso_energy_operator(self)

    def dense_per_material(self) -> dict[int, np.ndarray]:
        r"""Per-material operator matrix :math:`y\,\Sigma_{c,0}^{T}` (``[g_to, g_from]``)
        — the STORAGE-SIDE oracle view, not a production consumption mode.

        Returns ``{mid: M}`` with ``M @ φ_cell == apply(φ)_cell``, i.e.
        ``M = y · p0.T`` (``p0`` is stored ``[g_from, g_to]``), read
        DIRECTLY off the stored cross sections through the kernel's own
        :meth:`~orpheus.transport.kernels.TransferKernel.emission_matrix`
        — structurally independent of the
        :meth:`~orpheus.transport.material_field.TransferMaterialField.add_p0_source`
        einsum that realizes :meth:`apply`.  That independence is this
        method's job: the verification gates use it as the
        transpose-convention oracle for ``apply``/``apply_transpose``
        (`vv-principles` L11 — the two sides of an oracle pair must not
        share a realization).  Production materialization — the LHS-fold
        ``A = C − K_iso`` the homogeneous solver densifies — goes through
        the operator's own
        :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`
        apply-to-basis instead (taxonomy §12 step 5, which retired this
        docstring's earlier claim that the fold consumers read THIS
        method — they never did).  Each entry is a fresh copy.
        """
        return {
            mid: kernel.emission_matrix()
            for mid, kernel in self.transfer.per_material.items()
        }


class IsotropicScattering(IsotropicTransfer):
    r"""The P0 isotropic in-scatter energy operator :math:`\Sigma_{s,0}` — the :math:`S` term's energy binding.

    The scattering instance of :class:`IsotropicTransfer` (yield 1): a
    role class carrying only its extraction classmethod and its name.
    Consumed directly by the diffusion, homogeneous and CP energy paths
    and minted by
    :attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_energy`
    as the SN fast path's P0 half.
    """

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField", *, space: "FunctionSpace",
    ) -> "IsotropicScattering":
        r"""Tier-2 extract-and-mint: the P0 head of the facade's
        scattering channel, endomorphic on one ``space=``."""
        return cls(
            TransferMaterialField.scattering(mat_xs).at_order(0),
            domain=space,
            codomain=space,
        )


class IsotropicN2N(IsotropicTransfer):
    r"""The :math:`(n,2n)` isotropic energy operator :math:`2\,\Sigma_{2n}` — the :math:`N_{2n}` term's energy binding.

    The :math:`(n,2n)` instance of :class:`IsotropicTransfer` (yield
    :data:`~orpheus.transport.kernels.N2N_MULTIPLICITY`): a role class
    carrying only its extraction classmethod and its name. A DISTINCT
    *multiplication* term of the algebra (each event emits two neutrons;
    the yield also feeds the keff production accounting), summed with
    :class:`IsotropicScattering` for the isotropic emission source and
    minted by
    :attr:`~orpheus.transport.operators.n2n.N2NOperator.isotropic_energy`
    as the SN fast path's P0 half.

    **On the name.** The prefix names this class's ROLE — the
    :math:`\ell = 0` energy binding the scalar consumers (diffusion,
    homogeneous, CP, the SN k-balance and the ray seed's K_iso leaf)
    want — not a property of the reaction. `[M]` 2026-08-31 the GENDF
    files ORPHEUS ships store NL = 7 Legendre moments for MT=16, the
    same order as elastic; since #426 step 2 (2026-09-04) those moments
    reach the solve through the ANGULAR binding
    :class:`~orpheus.transport.operators.n2n.N2NOperator`, exactly as
    the elastic ones do through its scattering sibling. Until then this
    class's scalar domain was the whole model — a :math:`P_0`
    truncation the theory page records as history
    (``docs/theory/methods/sn/adjoint.rst`` §sn-n2n-p0-truncation).
    """

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField", *, space: "FunctionSpace",
    ) -> "IsotropicN2N":
        r"""Tier-2 extract-and-mint: the P0 head of the facade's
        :math:`(n,2n)` channel, endomorphic on one ``space=``."""
        return cls(
            TransferMaterialField.n2n(mat_xs).at_order(0),
            domain=space,
            codomain=space,
        )


@dataclass(eq=False)
class IsotropicFission(BoundOperator):
    r"""The fission energy operator :math:`F = |\chi\rangle\langle\nu\Sigma_f|` on the scalar flux.

    Per cell, :math:`(F\phi)_g = \chi_g \sum_{g'}\nu\Sigma_{f,g'}\,
    \phi_{g'}` — the rank-1 dyad: contraction against the production
    co-vector, emission along the spectrum. The ENERGY binding of the
    fission channel — what the k-eigenvalue outer iteration, the
    homogeneous :math:`K = A^{-1}F` composition, and the diffusion
    scalar composite consume. The ANGULAR binding of the same datum
    (the frame's :math:`\ell=0` conjugation of this dyad) is
    :class:`~orpheus.transport.operators.fission.FissionOperator`,
    which retains an instance of THIS class as its energy factor —
    the relation the transfer family spells as
    :attr:`~orpheus.transport.operators.transfer.TransferOperator.isotropic_energy`.

    **On the name.** Fission emission is isotropic BY CONSTRUCTION
    (:math:`\chi` carries no angular dependence — there is no
    anisotropic sibling), so the prefix carries family signal, not a
    contrast: it names this class's role in the greppable
    ``Isotropic*`` energy-binding family alongside
    :class:`IsotropicScattering` and :class:`IsotropicN2N` (CS4c
    step-4 design round §16.3).

    **Arithmetic home & caching.** Unlike its two loss-side siblings
    (which route per-material accumulation verbs), fission's arithmetic
    is the :class:`~orpheus.numerics.operator.RankOneOperator` dyad on
    CELLWISE factors, cached once via :attr:`kernel` (the CS4c step-3
    satellite ruling: composites are minted once per binding, not per
    apply — the read-through/depletion semantics of the pre-step-4
    ``FissionOperator.kernel`` are deliberately dropped; a depletion
    update re-binds the operator, as :math:`C` already does for
    :math:`\sigma_t`). The factors come from the retained
    :class:`~orpheus.transport.material_field.FissionMaterialField`
    (validated per-material :class:`~orpheus.transport.kernels.FissionKernel`
    pairs — the χ simplex/null law holds by construction), gathered
    over the binding's own bulk shape: SPACE FIRST, the ends size the
    data.

    **No assembly mode** — deliberate omission, not an oversight: the
    loss-side sibling assembles because the within-group stencil
    folds it into :math:`A`; fission is the eigenvalue OUTER source
    (within-group fission is zero) and no stencil consumer exists.
    Landing :meth:`assemble` here would be a gate-less capability
    (``coding-standards`` defer-until-consumer).

    Parameters
    ----------
    fission : FissionMaterialField
        The fission channel's kernel field — validated
        :math:`(\chi, \nu\Sigma_f)` pairs over the mesh layout.
    domain, codomain : FunctionSpace
        The two mandatory ends (kw-only, write-once): the scalar-flux
        space — plain ``(ng, *spatial)`` (SN k-outer, homogeneous) or
        the scalar-bulk composite ``FullFieldSpace`` (diffusion; the
        gather reads the interior). An ANGULAR composite is refused at
        construction — that consumer wants the angular binding,
        :class:`~orpheus.transport.operators.fission.FissionOperator`.
    """

    fission: "FissionMaterialField"
    # Fission emission is volumetric — bulk only, no face-trace action.
    # Class-level constant (unannotated: not a dataclass field).
    block_role = BlockRole.BULK

    def __post_init__(self) -> None:
        # CS4a K2 → CS4c per-END form (one shared guard, per end).
        self._assert_energy_extent_both_ends(
            self.fission.ng, operator="IsotropicFission",
        )
        # The scalar-bulk admission: the gather below sizes the factors
        # from the bulk shape, so a bulk whose LEADING axis is not the
        # group axis (an angular composite) must refuse here, loudly.
        bulk = self._bulk_scalar_space()
        if int(bulk.shape[0]) != int(self.fission.ng):
            raise TypeError(
                f"IsotropicFission requires a scalar (ng, *spatial) bulk "
                f"(leading axis = the group axis, ng="
                f"{self.fission.ng}); got bulk shape {tuple(bulk.shape)}. "
                f"An angular composite wants the ANGULAR binding — "
                f"FissionOperator."
            )

    def _bulk_scalar_space(self) -> "FunctionSpace":
        """The scalar (ng, *spatial) space the gathered factors ride —
        the composite's interior, or the plain scalar domain itself."""
        from orpheus.numerics.spaces.full_field_space import FullFieldSpace

        domain = self.domain
        if isinstance(domain, FullFieldSpace):
            interior = domain.interior_space
            if interior is None:
                raise TypeError(
                    "IsotropicFission: the bound composite domain carries "
                    "no interior space to size the gathered factors."
                )
            return interior
        return domain

    @classmethod
    def from_material_xs(
        cls, mat_xs: "MaterialXSField", *, space: "FunctionSpace",
    ) -> "IsotropicFission":
        r"""Tier-2 extract-and-mint: the facade's fission channel as
        validated per-material kernels, endomorphic on one ``space=``."""
        return cls(
            FissionMaterialField.from_material_xs(mat_xs),
            domain=space,
            codomain=space,
        )

    @property
    def data_ng(self) -> int:
        """The bound data's group count (family symmetry with the
        loss-side sibling's assembly read)."""
        return self.fission.ng

    @property
    def is_adjointable(self) -> bool:
        # apply_transpose realises the dual dyad |νΣf⟩⟨χ| (the factor
        # swap, a theorem of the rank-1 structure); is_invertible
        # inherits base False — a rank-1 production operator is singular.
        return True

    @cached_property
    def production_rate(self) -> "ReactionRateFunctional":
        r"""The production-rate co-vector :math:`\langle\nu\Sigma_f,\cdot\rangle`
        — fission's rank-1 ROW-FACTOR, as a typed inspectable
        :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
        (``coding-elegance`` Pattern 3: the per-cell fission source
        density is the most physically central diagnostic in
        criticality). The SAME type over :math:`\Sigma_a` is the
        absorption rate; ``k = production/absorption`` is their ratio.
        Cached once (see the class docstring's caching paragraph)."""
        from orpheus.transport.fields.cross_section_field import (
            CrossSectionField,
        )
        from orpheus.transport.reaction_rate_functional import (
            ReactionRateFunctional,
        )

        bulk = self._bulk_scalar_space()
        nu_sig_f = self.fission.gather_nu_sig_f(tuple(bulk.shape[1:]))
        return ReactionRateFunctional(
            CrossSectionField(values=nu_sig_f, space=bulk),
        )

    @cached_property
    def kernel(self) -> "TensorProductOperator":
        r"""The rank-1 TP kernel — the §5.6 integral structure and the ONE
        arithmetic home of every arm.

        ``outer(χ, production_rate) & IdentityOperator()``: the rank-1
        dyad on the gathered cellwise factors tensored with the
        spatial-axis broadcast, exactly the pre-step-4
        ``FissionOperator.kernel`` composition (the
        :meth:`~orpheus.numerics.operator.RankOneOperator.apply`
        reduction order is unchanged — bit-identity is gated). Its
        transpose IS the dual dyad :math:`|\nu\Sigma_f\rangle\langle\chi|`
        by the factor-swap theorem — no fission-specific transpose
        arithmetic exists anywhere."""
        from orpheus.numerics.operator import IdentityOperator, outer

        bulk = self._bulk_scalar_space()
        chi = self.fission.gather_chi(tuple(bulk.shape[1:]))
        return outer(chi, self.production_rate) & IdentityOperator()

    def apply(self, phi: "np.ndarray | FullField | object") -> "np.ndarray | FullField":
        r""":math:`\chi\,(\nu\Sigma_f\cdot\phi)` — the per-cell fission source (no :math:`1/k`).

        Bare ndarray / flux-carrier in → bare ndarray out (the
        model-portable contract). A scalar :class:`FullField` composite
        in → the closed scalar source composite out (the shared
        :func:`_scalar_composite_source`; the kernel is the SAME dyad
        either way). The :math:`1/k` division stays at the eigenvalue
        layer — this is a *linear* operator.
        """
        if isinstance(phi, FullField):
            return _scalar_composite_source(self, phi)
        return self.kernel.apply(_values_of(phi))

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r"""The dual dyad :math:`F^{T}\psi^* = \nu\Sigma_f\,(\chi\cdot\psi^*)`
        — the χ↔νΣf factor swap, routed through :attr:`kernel`'s
        :class:`~orpheus.numerics.operator.TensorProductOperator`
        transpose (single source; no fission-side arithmetic). This is
        the bare Euclidean :math:`F^{T}`; the metric Hilbert adjoint is
        the ``.H`` wrapper's job, composed from the spaces' own Riesz
        legs (CS4c step-4 §16.2 — nothing hand-rolled here).

        Bare-ndarray surface only, mirroring the family: the composite
        transpose belongs to the ANGULAR binding
        (:meth:`FissionOperator.apply_transpose`).
        """
        if isinstance(chi, FullField):
            raise TypeError(
                "IsotropicFission.apply_transpose: composite FullField "
                "transpose is not wired here — the composite adjoint "
                "belongs to the angular binding (FissionOperator); pass "
                "the bulk values."
            )
        return self.kernel.apply_transpose(_values_of(chi))

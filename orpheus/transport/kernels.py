r"""Interaction kernels — representation-free per-material physics data.

**The kernel/binding doctrine (Campaign 1, ruling 2 — "operators are born
bound").** An interaction *kernel* is the physics of ONE material's channel —
a Legendre transfer stack :math:`\{\Sigma_{c,\ell}\}` with its yield
:math:`y_c`, the fission factor pair :math:`(\chi, \nu\Sigma_f)` — with **no
representation attached**: no space, no mesh, no ``apply``. A *binding*
(a constructor: kernel × space [× assignment]) realizes the kernel as a
fully-bound operator with one domain, one codomain, one ``apply``. The
kernel is therefore the datum that survives re-realization: the same
``TransferKernel`` binds to an energy-only quotient space (the
homogeneous solver), to a discrete-ordinates composite (SN), or to a
moment space (a P\ :sub:`N`/moment solver) — polymorphism moves from
apply-time dispatch to construction time.

**One transfer family, the yield a datum (#426 step 2, ruled 2026-09-03).**
The tape stores, for EVERY channel, the same object — GENDF MF=6 holds
:math:`\sigma_{c,\ell}(g'\!\to\!g) = \sigma_c\,y_c\,f_{c,\ell}` (ENDF-102
Eq. 6.1/6.3, NJOY Eq. 242): a Legendre transfer stack carrying the
reaction's yield. Elastic has :math:`y = 1`, :math:`(n,2n)` has
:math:`y = 2`, and the two differ in **nothing else** — same record, same
convention, same non-separable :math:`g'\!\to\!g` structure. ORPHEUS keeps
the REACTION stack (the ingest divides the yield out, #427) and carries
:math:`y` beside it as :attr:`TransferKernel.multiplicity`: the loss side
counts :math:`\Sigma_c` once (one neutron absorbed per event —
``Mixture.absorption_xs``) while every emission verb scales by :math:`y`.
Scattering and :math:`(n,2n)` are two INSTANCES of this one kernel type,
not two types (`coding-standards` type-vs-property: isomorphic
realizations under one morphism — scale by :math:`y`); the production
accounting a :math:`y > 1` channel adds, :math:`(y-1)\,\Sigma_{c,0}^T\phi`,
vanishes for scattering by arithmetic, so no branch anywhere names the
channel. The channel constant :math:`y_{2n}` has ONE home,
:data:`N2N_MULTIPLICITY`.

**Why kernels view a** :class:`~orpheus.data.macro_xs.mixture.Mixture`.
The mixture is the honest provenance of every channel datum (its ``SigS``
/ ``Sig2`` sparse Legendre stacks and ``chi`` / ``SigP`` vectors are the
upstream nuclear-data reduction). The kernels copy that data **once, at
construction, into read-only arrays** — deliberately *not* aliasing any
carrier cache: measured 2026-08-20 (CS4a fact F4), the shipped
:meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.sig_s_legendre`
returned the production cache object itself, writable, so a consumer
mutation reached the loss matrix. A kernel datum is immutable by
construction (frozen dataclass + read-only buffers), which closes that
hazard for every kernel consumer — and since CS4a-R (EE-4) the producer
cache itself is frozen at build, closing it for the ~20 live carrier
consumers too; the kernel copy remains the guarantee that survives any
carrier rework.

**Storage conventions** (one home — this docstring; every array in this
module follows it):

* Matrix channels (:attr:`TransferKernel.moments`) are indexed
  ``[g_from, g_to]`` — the :attr:`Mixture.SigS` storage convention. The
  *action* on a flux is the transpose contraction :math:`(\Sigma^T\phi)_g
  = \sum_{g'} \Sigma(g'\!\to\!g)\,\phi_{g'}`; emission-side views
  (:meth:`TransferKernel.emission_matrix`, :meth:`FissionKernel.dyad`)
  return ``[g_to, g_from]`` operator matrices with ``M @ phi == action``.
* The Legendre index :math:`\ell` **is the Funk–Hecke eigenbasis index**:
  on the sphere, a rotationally-invariant transfer kernel is diagonal in
  spherical harmonics with eigenvalue :math:`\Sigma_{c,\ell}` on the whole
  degree-:math:`\ell` eigenspace. The stack is therefore the kernel's
  *spectral* representation, and CS2's frame-at-binding mints the angular
  frame from exactly this eigenbasis structure × the space's angular
  measure. The :math:`(2\ell+1)` of the addition theorem lives on the
  basis's reconstruction, never in a stored moment.
* The fission kernel stores the **factors** :math:`(\chi, \nu\Sigma_f)`,
  never the outer product: the dyad :math:`F = |\chi\rangle\langle
  \nu\Sigma_f|` is rank-1, its transpose is the factor swap *by theorem*
  (:math:`\mathrm{outer}(a,b)^T = \mathrm{outer}(b,a)`), and condensation
  couples the factors through the flux weighting — all three are
  statements about the factors that a stored matrix would obscure.

**What does NOT live here.** No mesh, no assignment (cell → material) —
that is the CS4b/CS1.5′ Medium layer's; no frames (CS2); no import from
the operator modules or the frame machinery (the C8 fence — the operator
tier points at the kernel tier, never the other way around).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import cached_property
from typing import TYPE_CHECKING, Final

import numpy as np

from orpheus.data.emission_spectrum import enforce_emission_spectrum

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture

__all__ = ["N2N_MULTIPLICITY", "TransferKernel", "FissionKernel"]

#: Neutrons emitted per :math:`(n,2n)` reaction — definitionally 2, a physics
#: constant of the channel, and its ONE home in the tree (XD-2): the solver
#: families read it from here, the transfer kernel of the channel carries it
#: as :attr:`TransferKernel.multiplicity`, and the data layer divides the
#: tape's yield out without ever naming it
#: (:func:`~orpheus.data.micro_xs.gendf._strip_transfer_yield`).
N2N_MULTIPLICITY: Final[int] = 2


def _read_only_copy(values: np.ndarray, *, name: str, ndim: int) -> np.ndarray:
    """A fresh, C-contiguous, float64, WRITE-PROTECTED copy of ``values``.

    The ``+ 0.0`` forces a new buffer even when ``ascontiguousarray``
    would pass the input through (the house idiom —
    ``EnergyAxis.__post_init__``), so freezing can never reach back and
    lock a caller's array, and no caller retains write access to kernel
    data.
    """
    arr = np.ascontiguousarray(values, dtype=float) + 0.0
    if arr.ndim != ndim:
        raise ValueError(
            f"{name} must be rank {ndim}; got shape {arr.shape}"
        )
    arr.setflags(write=False)
    return arr


@dataclass(frozen=True)
class TransferKernel:
    r"""The Legendre transfer stack :math:`\{\Sigma_{c,\ell}\}_{\ell=0}^{L}` of one material's channel, with its yield.

    ``moments[l]`` is the dense ``(ng, ng)`` group-transfer matrix of
    Legendre order :math:`\ell`, indexed ``[g_from, g_to]`` (module
    convention) — the REACTION stack :math:`\sigma_c f_{c,\ell}`, no
    yield folded in. :attr:`multiplicity` is the channel's yield
    :math:`y_c`, the number of neutrons one reaction emits: ``1`` for
    scattering (the identity — every scattering verb is unchanged by it)
    and :data:`N2N_MULTIPLICITY` for :math:`(n,2n)`. The stack length is
    the stored order plus one (:attr:`order` = :math:`L`), and the
    :math:`\ell`-index is the Funk–Hecke eigenbasis index (module
    docstring) — which is why *re-ordering* is the kernel's one
    representation-free morphism (:meth:`at_order`): dropping trailing
    eigenvalues, or reading zeros where the evaluation stored none, is
    meaningful before any space is in the room, while any rotation-,
    mesh- or ordinate-dependent operation is the *binding's*.

    Immutable by construction: the dataclass is frozen and every moment
    is a read-only fresh copy (``__post_init__``), so no consumer can
    reach production data through a kernel — the writable-cache-alias
    hazard this class exists to close. ``dataclasses.replace`` re-runs
    the validation (:meth:`at_order` routes through it); ``ng`` and
    ``order`` are derived properties, not fields, so
    ``replace(kernel, ng=...)`` raises ``TypeError`` by construction.
    """

    moments: tuple[np.ndarray, ...]
    multiplicity: int = 1

    def __post_init__(self) -> None:
        if len(self.moments) == 0:
            raise ValueError(
                "TransferKernel needs at least the l=0 moment; got an "
                "empty stack"
            )
        if (
            isinstance(self.multiplicity, bool)
            or not isinstance(self.multiplicity, int)
            or self.multiplicity < 1
        ):
            raise ValueError(
                f"TransferKernel.multiplicity is the number of neutrons one "
                f"reaction emits — a positive integer (1 for scattering, "
                f"{N2N_MULTIPLICITY} for (n,2n)); got {self.multiplicity!r}"
            )
        frozen = tuple(
            _read_only_copy(m, name=f"moments[{l}]", ndim=2)
            for l, m in enumerate(self.moments)
        )
        ng = frozen[0].shape[0]
        for l, m in enumerate(frozen):
            if m.shape != (ng, ng):
                raise ValueError(
                    f"every Legendre moment must be the same square "
                    f"(ng, ng) transfer matrix; moments[0] is "
                    f"{frozen[0].shape} while moments[{l}] is {m.shape}"
                )
        object.__setattr__(self, "moments", frozen)

    @classmethod
    def scattering(cls, mixture: "Mixture") -> "TransferKernel":
        """The scattering channel of ``mixture`` — a fresh copy of every
        sparse ``SigS[l]``, densified, yield 1; nothing is aliased."""
        return cls(
            moments=tuple(np.asarray(s.todense()) for s in mixture.SigS),
            multiplicity=1,  # elastic: one neutron out per event
        )

    @classmethod
    def n2n(cls, mixture: "Mixture") -> "TransferKernel":
        r"""The :math:`(n,2n)` channel of ``mixture`` — a fresh copy of
        every sparse ``Sig2[l]`` (the whole stack the tape stores, kept
        by the ingest since #426 step 1), densified, yield
        :data:`N2N_MULTIPLICITY`; nothing is aliased."""
        return cls(
            moments=tuple(np.asarray(s.todense()) for s in mixture.Sig2),
            multiplicity=N2N_MULTIPLICITY,
        )

    @property
    def ng(self) -> int:
        """Number of energy groups (derived from the stack; not a field)."""
        return self.moments[0].shape[0]

    @property
    def order(self) -> int:
        r"""The stored Legendre order :math:`L` (the stack holds :math:`L+1` moments)."""
        return len(self.moments) - 1

    @property
    def p0(self) -> np.ndarray:
        r"""The isotropic transfer matrix :math:`\Sigma_{c,0}` — ``moments[0]``, ``[g_from, g_to]``.

        The slice the isotropic energy binding consumes: the operator
        matrix of
        :meth:`~orpheus.transport.operators.isotropic_transfer.IsotropicTransfer.dense_per_material`
        is exactly :math:`y\,\mathtt{p0}^T` (:meth:`emission_matrix`).

        Aliasing semantics (CS4a-R EE-9c): this property returns the
        STORED read-only array itself (``p0 is moments[0]``) — unlike
        :meth:`emission_matrix` and :meth:`FissionKernel.dyad`, which
        mint a fresh WRITABLE copy per call. Do not cache one expecting
        the other's semantics.
        """
        return self.moments[0]

    @cached_property
    def is_isotropic(self) -> bool:
        r"""``True`` iff every moment above :math:`\ell = 0` is exactly zero —
        the channel's :math:`\Lambda_{\ell\ge 1}` is the zero operator.

        A structural fact of the datum, not a tolerance: an absent
        section, an ``NL = 1`` evaluation, or a stack padded by
        :meth:`at_order` are exactly zero above :math:`\ell = 0`, and the
        angular binding uses this to skip the :math:`R\Lambda M` product
        that would return exact zeros. Cached (the stack is immutable).
        """
        return not any(m.any() for m in self.moments[1:])

    def emission_matrix(self) -> np.ndarray:
        r"""The emission-side operator matrix :math:`y\,\Sigma_{c,0}^{T}` — ``[g_to, g_from]``.

        ``M @ phi`` is the per-cell isotropic emission of the channel:
        the P0 in-scatter for :math:`y = 1`, the :math:`(n,2n)` source for
        :math:`y = 2`. A fresh copy per call (the storage-side-view
        convention of
        :meth:`~orpheus.transport.operators.isotropic_transfer.IsotropicTransfer.dense_per_material`,
        which this equals entry for entry). The multiplicity enters the
        emission HERE and in the field's verbs — never the stored stack.
        """
        return np.ascontiguousarray(self.multiplicity * self.p0.T)

    def at_order(self, order: int) -> "TransferKernel":
        r"""The P\ :sub:`order` view of the stack — the first ``order + 1``
        moments, with **zero moments above the stored order**.

        The identity at ``order == self.order`` (returns ``self``; the
        kernel is immutable). Below it, the trailing eigenvalues are
        dropped — the binding's truncation. Above it, the missing moments
        are exactly zero, which is the EVALUATION's statement, not an
        invention: a stack shorter than the request is complete — an
        absent section, or an ``NL = 1`` section declaring isotropy — and
        a consumer that reads its :math:`\ell \ge 1` moments reads the
        zeros the tape means. ⚠ The one stack that must never be padded —
        a channel at GROUPR's cap, asked for a moment the evaluation
        capped away — is a CALLER's obligation, not this method's: the
        kernel cannot tell a complete short stack from a capped one (the
        tape gives no such datum), and the caller that keeps the promise
        today is the SN solve, whose clamp reads the SCATTERING stack
        alone (ruling O-1), so a request above a shipped stack's order
        only ever reaches the shorter sibling channel (plan §4.3, ruled
        2026-09-03). A second caller that pads owes the same argument.
        """
        if order < 0:
            raise ValueError(
                f"TransferKernel.at_order requires order >= 0; got {order}"
            )
        if order == self.order:
            return self
        if order < self.order:
            return replace(self, moments=self.moments[: order + 1])
        zero = np.zeros((self.ng, self.ng))
        return replace(
            self,
            moments=(*self.moments, *((zero,) * (order - self.order))),
        )


@dataclass(frozen=True)
class FissionKernel:
    r"""The fission channel of one material — the factor pair :math:`(\chi, \nu\Sigma_f)`.

    The fission source is the rank-1 dyad :math:`F = |\chi\rangle\langle
    \nu\Sigma_f|`: contraction of the flux against the production
    co-vector :math:`\nu\Sigma_f`, then emission along the spectrum
    :math:`\chi`. The kernel stores the **factors, never the outer
    product** — rank-1-ness, the transpose-is-the-factor-swap theorem
    (gated as algebra: ``dyad().T == outer(nu_sig_f, chi)``), and the
    χ↔νΣf-coupled condensation are all facts *about the factors* that a
    stored ``(ng, ng)`` matrix would obscure. On an energy-only space the
    dyad fully determines the adjoint (counting measure ⟹ metric adjoint
    ≡ transpose ≡ the swap); on angular/moment spaces the *binding*
    supplies the retraction/embedding pair from the space's measure (the
    F4-addendum composite) — neither is this class's concern, which is
    exactly the point.

    The χ simplex/null law is enforced through the one shared home
    (:func:`~orpheus.data.emission_spectrum.enforce_emission_spectrum`,
    the same law ``Mixture`` and ``Isotope`` run): a producing kernel
    (:math:`\nu\Sigma_f > 0` anywhere) must carry a probability simplex,
    a non-producing one the null spectrum. ``dataclasses.replace``
    re-runs it (G1.8: replace re-validates).

    **Consumer status (Q2; updated at CS4c step 4 — the rebind this
    paragraph awaited LANDED).** The production fission path consumes
    this datum: every material's pair enters a
    :class:`~orpheus.transport.material_field.FissionMaterialField`
    (validated per material by THIS constructor), whose gathered factors
    feed the energy binding
    :class:`~orpheus.transport.operators.isotropic_transfer.IsotropicFission`
    (the k-outer / homogeneous / diffusion realization) and, through it,
    the angular composite
    :class:`~orpheus.transport.operators.fission.FissionOperator` (the
    frame's :math:`\ell=0` conjugation). The χ↔νΣf-coupled condensation
    is gated (G-F1, ``tests/transport/test_kernels.py``).

    Fission stays its own kernel type under the type-vs-property rule the
    transfer family was ruled by: its realization is a separable rank-1
    dyad (not isomorphic to a Legendre stack) under a different morphism
    (the eigenvalue's :math:`1/k`).
    """

    chi: np.ndarray
    nu_sig_f: np.ndarray

    def __post_init__(self) -> None:
        chi = _read_only_copy(self.chi, name="chi", ndim=1)
        nu_sig_f = _read_only_copy(self.nu_sig_f, name="nu_sig_f", ndim=1)
        if chi.shape != nu_sig_f.shape:
            raise ValueError(
                f"the fission factors are two (ng,) vectors of one "
                f"material; got chi {chi.shape} vs nu_sig_f "
                f"{nu_sig_f.shape}"
            )
        chi = np.asarray(
            enforce_emission_spectrum(
                chi, is_producing=bool(np.any(nu_sig_f > 0))
            ),
            dtype=float,
        )
        chi.setflags(write=False)
        object.__setattr__(self, "chi", chi)
        object.__setattr__(self, "nu_sig_f", nu_sig_f)

    @classmethod
    def from_mixture(cls, mixture: "Mixture") -> "FissionKernel":
        r"""The kernel of ``mixture``'s fission channel: :math:`(\chi,\ \nu\Sigma_f = \mathtt{SigP})`."""
        return cls(chi=mixture.chi, nu_sig_f=mixture.SigP)

    @property
    def ng(self) -> int:
        """Number of energy groups (derived from the factors; not a field)."""
        return self.chi.shape[0]

    def dyad(self) -> np.ndarray:
        r"""The dense dyad :math:`|\chi\rangle\langle\nu\Sigma_f|` — ``[g_to, g_from]``, fresh per call.

        ``dyad() @ phi`` is the fission source. A *view* for oracles and
        small dense consumers — production realizations keep the factors
        (rank-1 apply is :math:`\chi\,(\nu\Sigma_f\cdot\phi)`, never a
        matmul against this).
        """
        return np.outer(self.chi, self.nu_sig_f)

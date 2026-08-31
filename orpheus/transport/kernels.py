r"""Interaction kernels — representation-free per-material physics data.

**The kernel/binding doctrine (Campaign 1, ruling 2 — "operators are born
bound").** An interaction *kernel* is the physics of ONE material's channel —
the Legendre transfer stack :math:`\{\Sigma_{s\ell}\}`, the :math:`(n,2n)`
matrix, the fission factor pair :math:`(\chi, \nu\Sigma_f)` — with **no
representation attached**: no space, no mesh, no ``apply``. A *binding*
(a constructor: kernel × space [× assignment]) realizes the kernel as a
fully-bound operator with one domain, one codomain, one ``apply``. The
kernel is therefore the datum that survives re-realization: the same
``ScatteringKernel`` binds to an energy-only quotient space (the
homogeneous solver), to a discrete-ordinates composite (SN), or to a
moment space (a P\ :sub:`N`/moment solver) — polymorphism moves from
apply-time dispatch to construction time.

**Why kernels view a** :class:`~orpheus.data.macro_xs.mixture.Mixture`.
The mixture is the honest provenance of every channel datum (its ``SigS``
/ ``Sig2`` sparse matrices and ``chi`` / ``SigP`` vectors are the
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

* Matrix channels (:attr:`ScatteringKernel.moments`,
  :attr:`N2NKernel.matrix`) are indexed ``[g_from, g_to]`` — the
  :attr:`Mixture.SigS` storage convention. The *action* on a flux is the
  transpose contraction :math:`(\Sigma^T\phi)_g = \sum_{g'}
  \Sigma(g'\!\to\!g)\,\phi_{g'}`; emission-side views
  (:meth:`N2NKernel.emission_matrix`, :meth:`FissionKernel.dyad`) return
  ``[g_to, g_from]`` operator matrices with ``M @ phi == action``.
* The Legendre index :math:`\ell` **is the Funk–Hecke eigenbasis index**:
  on the sphere, a rotationally-invariant transfer kernel is diagonal in
  spherical harmonics with eigenvalue :math:`\Sigma_{s\ell}` on the whole
  degree-:math:`\ell` eigenspace. The stack is therefore the kernel's
  *spectral* representation, and CS2's frame-at-binding mints the angular
  frame from exactly this eigenbasis structure × the space's angular
  measure.
* The fission kernel stores the **factors** :math:`(\chi, \nu\Sigma_f)`,
  never the outer product: the dyad :math:`F = |\chi\rangle\langle
  \nu\Sigma_f|` is rank-1, its transpose is the factor swap *by theorem*
  (:math:`\mathrm{outer}(a,b)^T = \mathrm{outer}(b,a)`), and condensation
  couples the factors through the flux weighting — all three are
  statements about the factors that a stored matrix would obscure.

**What does NOT live here.** No mesh, no assignment (cell → material) —
that is the CS4b/CS1.5′ Medium layer's; no frames (CS2); no import from
``scattering.py`` or the frame machinery (the C8 fence —
``ScatteringOperator`` re-points to :class:`ScatteringKernel` at CS4c,
not the other way around). The :math:`(n,2n)` channel keeps its own
kernel rather than folding into scattering: it enters the **loss** side
asymmetrically (one neutron absorbed, two emitted — ``absorption_xs``
counts :math:`\Sigma_{2n}` once while emission carries the factor
:attr:`N2NKernel.multiplicity`), a twice-ruled physics decision
(``isotropic_scattering.py`` header).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.data.emission_spectrum import enforce_emission_spectrum

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture

__all__ = ["ScatteringKernel", "N2NKernel", "FissionKernel"]


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
class ScatteringKernel:
    r"""The Legendre transfer stack :math:`\{\Sigma_{s\ell}\}_{\ell=0}^{L}` of one material.

    ``moments[l]`` is the dense ``(ng, ng)`` group-transfer matrix of
    Legendre order :math:`\ell`, indexed ``[g_from, g_to]`` (module
    convention). The stack length is the truncation order plus one
    (:attr:`order` = :math:`L` of the P\ :sub:`L` expansion), and the
    :math:`\ell`-index is the Funk–Hecke eigenbasis index (module
    docstring) — which is why *truncation* is the kernel's one
    representation-free morphism: dropping trailing eigenvalues is
    meaningful before any space is in the room, while any rotation-,
    mesh- or ordinate-dependent operation is the *binding's*.

    Immutable by construction: the dataclass is frozen and every moment
    is a read-only fresh copy (``__post_init__``), so no consumer can
    reach production data through a kernel — the writable-cache-alias
    hazard this class exists to close. ``dataclasses.replace`` re-runs
    the validation (:meth:`truncated` routes through it); ``ng`` and
    ``order`` are derived properties, not fields, so
    ``replace(kernel, ng=...)`` raises ``TypeError`` by construction.
    """

    moments: tuple[np.ndarray, ...]

    def __post_init__(self) -> None:
        if len(self.moments) == 0:
            raise ValueError(
                "ScatteringKernel needs at least the l=0 moment; got an "
                "empty stack"
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
    def from_mixture(cls, mixture: "Mixture") -> "ScatteringKernel":
        """The kernel of ``mixture``'s scattering channel — a fresh copy of
        every sparse ``SigS[l]``, densified; nothing is aliased."""
        return cls(
            moments=tuple(np.asarray(s.todense()) for s in mixture.SigS)
        )

    @property
    def ng(self) -> int:
        """Number of energy groups (derived from the stack; not a field)."""
        return self.moments[0].shape[0]

    @property
    def order(self) -> int:
        r"""The Legendre truncation order :math:`L` (the stack holds :math:`L+1` moments)."""
        return len(self.moments) - 1

    @property
    def p0(self) -> np.ndarray:
        r"""The isotropic transfer matrix :math:`\Sigma_{s0}` — ``moments[0]``, ``[g_from, g_to]``.

        The slice the isotropic in-scatter pair consumes: the operator
        matrix of
        :meth:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering.dense_per_material`
        is exactly ``p0.T``.

        Aliasing semantics (CS4a-R EE-9c): this property returns the
        STORED read-only array itself (``p0 is moments[0]``) — unlike
        :meth:`N2NKernel.emission_matrix` and :meth:`FissionKernel.dyad`,
        which mint a fresh WRITABLE copy per call. Do not cache one
        expecting the other's semantics.
        """
        return self.moments[0]

    def truncated(self, order: int) -> "ScatteringKernel":
        r"""The P\ :sub:`order` sub-kernel — the first ``order + 1`` moments.

        The identity at ``order == self.order``. **Refuses**
        ``order > self.order``: the data carries no moments beyond
        :math:`L`, and inventing zero moments would silently misreport
        the material's anisotropy content as a measured zero (a fabricated
        datum — the campaign's O1 tell).
        """
        if not 0 <= order <= self.order:
            raise ValueError(
                f"ScatteringKernel.truncated requires 0 <= order <= "
                f"{self.order} (the data's own truncation); got {order}. "
                f"Moments beyond the stored order do not exist and are "
                f"not invented."
            )
        return replace(self, moments=self.moments[: order + 1])


@dataclass(frozen=True)
class N2NKernel:
    r"""The :math:`(n,2n)` channel of one material: :math:`\Sigma_{2n}` + the multiplicity.

    :attr:`matrix` is the dense ``(ng, ng)`` *reaction* matrix
    :math:`\Sigma_{2n}`, indexed ``[g_from, g_to]`` (module convention) —
    the raw cross section, with **no** multiplicity folded in. The factor
    lives once *in this module*, in :attr:`multiplicity` (`[M]` CS4a-R
    XD-2: the wider tree still spells the constant as a literal ``2`` at
    12 production sites across the solver families — their consolidation
    onto this one home is CS4c's rebind obligation), and the
    emission-side operator matrix :math:`2\,\Sigma_{2n}^{T}` is minted
    by :meth:`emission_matrix`. Keeping reaction and emission distinct is
    the loss-side channel ruling made structural: absorption counts
    :math:`\Sigma_{2n}` once (one neutron absorbed) while emission
    carries the 2 (two emitted) — a kernel that stored
    :math:`2\Sigma_{2n}` could not serve the loss side without dividing
    the physics constant back out.
    """

    matrix: np.ndarray

    #: Neutrons emitted per (n,2n) reaction — definitionally 2, a physics
    #: constant of the channel (a ClassVar, not a field: there is no
    #: legal (n,2n) kernel with another value, so it is unrepresentable).
    multiplicity: ClassVar[int] = 2

    def __post_init__(self) -> None:
        frozen = _read_only_copy(self.matrix, name="matrix", ndim=2)
        if frozen.shape[0] != frozen.shape[1]:
            raise ValueError(
                f"the (n,2n) matrix is a square (ng, ng) group-transfer "
                f"matrix; got shape {frozen.shape}"
            )
        object.__setattr__(self, "matrix", frozen)

    @classmethod
    def from_mixture(cls, mixture: "Mixture") -> "N2NKernel":
        """The kernel of ``mixture``'s :math:`(n,2n)` channel — a fresh
        densified copy of the sparse ``Sig2``; nothing is aliased."""
        return cls(matrix=np.asarray(mixture.Sig2.todense()))

    @property
    def ng(self) -> int:
        """Number of energy groups (derived from the matrix; not a field)."""
        return self.matrix.shape[0]

    def emission_matrix(self) -> np.ndarray:
        r"""The emission-side operator matrix :math:`2\,\Sigma_{2n}^{T}` — ``[g_to, g_from]``.

        ``M @ phi`` is the per-cell (n,2n) source. A fresh copy per call
        (the storage-side-view convention of
        :meth:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N.dense_per_material`,
        which this equals entry for entry). The multiplicity enters HERE
        and only here — :attr:`matrix` stays the raw reaction XS.
        """
        return np.ascontiguousarray(self.multiplicity * self.matrix.T)


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
    :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicFission`
    (the k-outer / homogeneous / diffusion realization) and, through it,
    the angular composite
    :class:`~orpheus.transport.operators.fission.FissionOperator` (the
    frame's :math:`\ell=0` conjugation). The χ↔νΣf-coupled condensation
    is gated (G-F1, ``tests/transport/test_kernels.py``).
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

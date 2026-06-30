r"""Multigroup scattering source as a :class:`LinearOperator`.

This module owns the **scattering source operator** :math:`S` from the
operator-algebra view of the Boltzmann transport equation,

.. math::

    (L - S - F)\,\psi = q

where :math:`L` is the streaming + collision operator, :math:`F` is the
fission operator, and :math:`S` here aggregates **all secondary-emission
channels that depend on the in-cell flux**:

* **P0 isotropic in-scatter** :math:`\Sigma_s^0(g'\to g)\,\phi_{g'}`,
* **Pℓ Galerkin reconstruction on real spherical harmonics**
  :math:`Y_\ell^m` for :math:`\ell\ge 1` from the angular flux moments,
* **(n,2n) doubling** :math:`2\,\Sigma_{2n}(g'\to g)\,\phi_{g'}`, which
  ORPHEUS folds into the scattering side because its cell-by-cell
  bookkeeping is identical to in-scatter (vectorise-by-material,
  add-into-:math:`Q`).

The operator advertises only :pydata:`CAP_APPLY`. There is no useful
:meth:`solve`: the discrete :math:`S` is rank :math:`O(N_{\text{cells}}
\cdot N_{\text{groups}})` with no tractable inverse — it is *applied*,
never *inverted*. This is the canonical example the capability-set
design (Wave A Issue 1) was built for: forcing :math:`S` to provide a
``solve`` stub that raises ``NotImplementedError`` is harmful, so we
declare what :math:`S` **can** do and let composers refuse to build
:math:`S^{-1}` at construction time.

Per Cardinal Rule 2 (architecture) this lifts the
``SNSolver._add_scattering_source`` / ``SNSolver._build_aniso_scattering``
/ ``SNSolver._add_n2n_source`` math out of the solver and into a
single operator object. The math is **moved verbatim** (Wave D Issue 13
is a bit-identical extraction); the only change is the home of the code.
:class:`SNSolver` retains thin delegator methods so the EigenvalueSolver
Protocol surface and the ``test_solver_components.py`` probe-tests
keep working unchanged.

Mathematical content
====================

P0 in-scatter
-------------

Per material region :math:`m`, with :math:`\Sigma_{s,0}^m \in
\mathbb{R}^{G\times G}` the row-from / column-to scattering matrix at
order :math:`\ell = 0`,

.. math::

    Q^{(0)}_g(\vec r)
    = \sum_{g'} \Sigma_{s,0}^m(g'\to g)\,\phi_{g'}(\vec r)
    = \bigl(\phi(\vec r)\,\cdot\,\Sigma_{s,0}^m\bigr)_g

with :math:`\phi(\vec r) \in \mathbb{R}^G` and the dot product in
matrix-row order. This is isotropic: it appears identically in every
ordinate via the :math:`1/W = 1/\sum_n w_n` normalisation that the
sweep applies.

Pℓ Galerkin projection on :math:`Y_\ell^m`
------------------------------------------

For anisotropic scattering of order :math:`L \ge 1`, the angular flux
is projected onto the real spherical harmonics
:math:`Y_\ell^m(\hat\Omega)`,

.. math::

    \phi^{\ell m}_g(\vec r)
    = \int_{4\pi}\psi_g(\vec r,\hat\Omega)\,Y_\ell^m(\hat\Omega)\,
      d\hat\Omega
    \;\approx\; \sum_n w_n\,\psi_{n,g}(\vec r)\,Y_\ell^m(\hat\Omega_n)

(the right-hand side is the discrete quadrature; canonical equation
label :eq:`flux-moments` lives in :ref:`theory-discrete-ordinates`).
The :math:`\ell\ge 1` contribution to the per-ordinate scattering
source uses the addition theorem
:math:`\sum_m Y_\ell^m(\hat\Omega)\,Y_\ell^m(\hat\Omega') =
P_\ell(\hat\Omega\cdot\hat\Omega')`,

.. math::

    Q^{(\ell\ge 1)}_{n,g}(\vec r)
    = \sum_{\ell=1}^{L}\,(2\ell+1)\,\sum_m Y_\ell^m(\hat\Omega_n)\,
      \sum_{g'}\Sigma_{s,\ell}^m(g'\to g)\,\phi^{\ell m}_{g'}(\vec r).

The :math:`(2\ell+1)` factor is the discrete-orthogonality
normalisation
:math:`\langle Y_\ell^m | Y_{\ell'}^{m'}\rangle = (4\pi/(2\ell+1))\,
\delta_{\ell\ell'}\delta_{mm'}` working out across both projection
and reconstruction. The Galerkin frame is real spherical harmonics, not
complex — see :class:`orpheus.sn.quadrature.LebedevSphere`.

(n,2n) doubling
---------------

The (n,2n) reaction emits **two** neutrons per absorption, producing a
secondary source :math:`2\,\Sigma_{2n}^m(g'\to g)\,\phi_{g'}`. Because
this depends on the scalar flux only and adds isotropically (in this
discretisation), it folds into the scattering side of the algebra:
:meth:`apply` returns the sum of P0 + (n,2n) + Pℓ contributions.

Capability advertisement
========================

:pydata:`capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})`. No
``solve`` (rank-deficiency on the :math:`\ell=0` block prevents efficient
inversion). The adjoint :math:`S^{T}` IS advertised (campaign #276 A2b /
`#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_):
:meth:`~ScatteringOperator.apply_transpose` rides the harmonic-frame
:attr:`~ScatteringOperator.full_scatter_kernel`, whose transpose
:math:`(R\circ(\Lambda+N_{2n})\circ M)^{T} = M^{T}\circ(\Lambda+N_{2n})^{T}
\circ R^{T}` falls out of
:meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose` for free
— so the WHOLE source (iso :math:`\ell=0` + aniso :math:`\ell\ge1` +
(n,2n)) transposes in one expression.

The §5.6 Kernel reading — scattering as an integral kernel
==========================================================

In the grand-report §5.6 suffix law (see
:mod:`orpheus.transport.operators.integral_kernel_operator`) scattering is a
**Kernel**: a *nonlocal* operator whose anisotropic action integrates
the flux against a measure on the angular axis (the :math:`P_\ell`
source at ordinate :math:`\hat\Omega_n` reads the flux at *every*
ordinate, via projection-then-reconstruction). #257 S6 makes this a
type-level fact — :class:`ScatteringOperator` exposes the integral
structure through a :attr:`~ScatteringOperator.kernel` property, the
typed :class:`~orpheus.numerics.operator.OperatorProduct`

.. math::

    \mathrm{kernel} \;=\; R \;\circ\; \Lambda_{\ell\ge 1} \;\circ\; M ,

so the operator satisfies the
:class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
Protocol. The :attr:`~ScatteringOperator.kernel` exposes the
anisotropic Legendre redistribution (the genuinely-nonlocal-in-angle
part); the isotropic :math:`\ell=0` :math:`P_0` in-scatter fast path
and the :math:`(n,2n)` doubling are the LOCAL / separate components of
the full :meth:`~ScatteringOperator.apply`. The ``kernel`` is the §5.6
*semantic reading* of the existing :math:`R\Lambda M` aniso path
(:meth:`~ScatteringOperator.build_aniso_source` /
:meth:`~ScatteringOperator._aniso_source_from_moment_values`); the 5
:meth:`~ScatteringOperator.apply` dispatch arms are UNCHANGED in S6
(additive, bit-identical), and the two compose the same numpy chain
(pinned at 0 ULP by
``tests/sn/operators/test_scattering_kernel_crosscheck.py``). The
producer-side :math:`1/W` normalisation (lesson L18) lives OUTSIDE the
kernel, at the :meth:`~ScatteringOperator.apply` boundary.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, singledispatchmethod
from typing import TYPE_CHECKING, Any, cast, overload

import numpy as np

from orpheus.numerics.frame import GalerkinFrame
from orpheus.transport.frames import HarmonicFrame
from orpheus.numerics.operator import (
    BlockRole,
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperator,
    OperatorProduct,
)

# Runtime imports of the flux / source types — required at module load
# time because :func:`singledispatchmethod.register` dispatches on the
# runtime class.  These three modules form a leaf in the SN package
# dependency graph (they do not import scattering.py), so the imports
# are circular-import-safe.
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import (
    ScalarSourceSink,
    AngularSourceSink,
    BoundarySourceSink,
    HarmonicMomentSourceSink,
)
from orpheus.transport.full_field import FullField
from orpheus.transport.timed_full_field import TimedFullField
# Runtime (not TYPE_CHECKING-only) — the windowed moment-iterate ``apply``
# arm registers on this type via ``@apply.register``, which needs the
# concrete class at class-definition time.  Circular-import-safe (the
# transport field types do not import scattering.py; same as ScalarFlux /
# AngularFlux above).
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces import FullFieldSpace


__all__ = ["LegendreMomentScattering", "ScatteringOperator"]


def _spatial_moments_of(phi: "np.ndarray | ScalarFlux") -> int:
    r"""Within-cell spatial-moment count per axis carried by ``phi``, or ``1``.

    The scattering-source accumulators (:meth:`ScatteringOperator._assemble_per_ordinate_source`)
    must match the driving flux's spatial-moment width so the in-place
    ``Σ_s ⊗ I`` accumulation does not broadcast-fail (#240 D5b-S3).  A typed
    :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` carries the width
    as a :class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
    factor on its space (the single source of truth — read it OFF the field, do
    not re-thread an int); a bare ``ndarray`` φ has no factor → ``1``."""
    return getattr(phi, "spatial_moments_per_axis", 1)


@dataclass(frozen=True)
class LegendreMomentScattering(LinearOperator):
    r"""Per-ℓ block-diagonal scattering on harmonic-moment space.

    The §15.2 / §10 sum-of-tensor-products form

    .. math::

        \Lambda \;=\; \sum_{\ell=1}^{L}\, \mathbf{P}_\ell \otimes \Sigma_{s,\ell}

    where :math:`\mathbf{P}_\ell` selects the :math:`\ell`-th harmonic
    block and :math:`\Sigma_{s,\ell}` is the per-material per-ℓ
    Legendre cross-section matmul on the energy-group axis.

    For an input moment field
    :math:`\phi_\ell^m(\vec r) \in \mathbb{R}^{(L+1) \times (2L+1) \times N_g \times N_x \times N_y}`
    (Issue #196 PR-INDEX-4 — principled ``g`` leading-after-moments),
    the action is

    .. math::

        (\Lambda \phi)_\ell^m(\vec r)\bigg|_{g}
        \;=\; \sum_{g'} \Sigma_{s,\ell}^{m(\vec r)}(g' \to g)\,
              \phi_\ell^m(\vec r)\bigg|_{g'},

    where :math:`m(\vec r)` is the material id at cell :math:`\vec r`
    (the per-material structure is folded into the cell axis via
    ``cells_by_mat``).

    The :math:`\ell = 0` block is **skipped by default** because the
    project's P0 isotropic in-scatter is handled by the separate
    :meth:`ScatteringOperator.add_iso_source` fast path (a reaction-rate
    formulation that does not go through the moment tensor). Setting
    ``skip_l0 = False`` includes the :math:`\ell = 0` block (used by
    the LinearOperator-surface :meth:`ScatteringOperator.apply` when a
    full :math:`R \Lambda M \psi` composition is required).

    Implementation
    --------------

    The two outer Python loops are over **materials** (1-10 typically)
    and **Legendre order ℓ** (0-3 typically) — both small structural
    dimensions, NOT smoking guns. Inside each ``(mid, ℓ)`` pair, ONE
    :func:`numpy.einsum` handles all :math:`m`, cells, and groups in
    a single contraction. The total flop count is identical to the
    legacy hand-rolled version inside
    :meth:`ScatteringOperator.build_aniso_source`; the iteration over
    ordinates that *was* there is no longer needed because moment-space
    scattering does not couple ordinates.

    Capability set: :pydata:`{"apply", "apply_transpose"}`. There is no
    efficient :meth:`solve` (the moment-space scattering matrix is
    rank-deficient on the :math:`\ell = 0` block by design). :meth:`apply_transpose`
    IS advertised (campaign #276): :math:`\Lambda^{T}` is the per-ℓ group-axis
    transpose of the block-diagonal :math:`\Sigma_{s,\ell}` matmul — the ONLY
    group-asymmetric factor of the frame-conjugated scattering kernel
    :math:`R\circ\Lambda\circ M`, so :math:`(R\circ\Lambda\circ M)^{T} =
    M^{T}\circ\Lambda^{T}\circ R^{T}` falls out of
    :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose` for free
    (the :math:`M`/:math:`R` face transposes landed in the Frame carve).

    Parameters
    ----------
    mat_xs : MaterialXSField
        Macroscopic XS field carrying both the per-material Legendre
        scattering matrices and the cell-to-material map.  Issue #197
        PR-TYPED-1 — the per-material dispatch loop lives ONLY inside
        :meth:`MaterialXSField.apply_legendre_scattering_moments`.
    L : int
        Maximum Legendre order :math:`L` retained.
    skip_l0 : bool, default ``True``
        Skip the :math:`\ell = 0` block (handled separately by P0
        in-scatter). Set to ``False`` for the full LinearOperator
        composition :math:`R \Lambda M \psi`.
    """

    mat_xs: "MaterialXSField"
    L: int
    skip_l0: bool = True
    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )

    @property
    def is_adjointable(self) -> bool:
        # Λ exposes its group-transpose Σ_{s,ℓ}^T (apply_transpose), so the
        # metric-aware .H is free; caps ⊇ CAP_APPLY_TRANSPOSE. is_invertible
        # inherits base False — a per-ℓ source map is not invertible.
        return True

    def apply(
        self, moments: "np.ndarray | HarmonicMomentFlux",
    ) -> "np.ndarray | HarmonicMomentSourceSink":
        r"""Apply :math:`\Lambda` to a moment field — the **role-changing** edge.

        :math:`\Lambda` (the per-:math:`\ell` group-transfer cross sections
        :math:`\Sigma_{s,\ell}`) maps a flux moment to the in-scatter **source**
        moment it emits: it is the sole role-changing edge of the
        ``(angular ⊗ moment) × (flux ⊗ source)`` carrier grid (axis-preserving,
        role flux → source). The typed arm makes that role change EXPLICIT in
        the signature.

        Parameters
        ----------
        moments : np.ndarray or HarmonicMomentFlux
            Flux moment field of shape ``(L+1, 2L+1, ng, *spatial)`` (Issue
            #196 PR-INDEX-4 — principled).  The :math:`m`-axis is the
            addition-theorem-shifted index where slot ``l + m`` holds
            the :math:`(\ell, m)` entry; entries outside
            :math:`|m| \le \ell` are conventionally zero.

            Typed
            :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
            is accepted; when supplied, the return is a typed
            :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
            (the scattered moment SOURCE) with matching ``L`` / ``mesh`` /
            spatial-moment width.  Bare ndarray is still accepted (the
            endomorphic moment-space view the :math:`R\circ\Lambda\circ M`
            :attr:`~ScatteringOperator.kernel` ``OperatorProduct`` composes on);
            its return is ndarray (same coefficient space, role implicit).

        Returns
        -------
        np.ndarray or HarmonicMomentSourceSink
            The scattered moment source, same shape as ``moments``.  The
            :math:`\ell = 0` block is zero when ``skip_l0`` is ``True``;
            otherwise the P0 in-scatter contribution is included.  Typed in →
            typed source out; ndarray in → ndarray out.

        Notes
        -----
        Per-material dispatch is encapsulated by
        :meth:`MaterialXSField.apply_legendre_scattering_moments` —
        Issue #197 PR-TYPED-1 collapses the ``for mid, (ix, iy) in
        cells_by_mat.items()`` loop into a typed verb. Both the typed and
        ndarray arms route through that SINGLE kernel (Pattern 2); they differ
        only in the carrier wrap.
        """
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
        if isinstance(moments, HarmonicMomentFlux):
            out_values = self.mat_xs.apply_legendre_scattering_moments(
                moments.values, L=self.L, skip_l0=self.skip_l0,
            )
            # flux moment → source moment: the explicit role change.
            return HarmonicMomentSourceSink.from_mesh_and_L(
                out_values, moments.mesh, moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.mat_xs.apply_legendre_scattering_moments(
            moments, L=self.L, skip_l0=self.skip_l0,
        )

    def apply_transpose(
        self, moments: "np.ndarray | HarmonicMomentSourceSink",
    ) -> "np.ndarray | HarmonicMomentFlux":
        r"""Apply :math:`\Lambda^{T}` — the per-ℓ group-transpose (the role-REVERSING edge).

        The bare Euclidean transpose of :meth:`apply`: where :math:`\Lambda` maps a
        flux moment to the in-scatter SOURCE moment it emits (flux → source),
        :math:`\Lambda^{T}` maps a source moment back into the flux-moment space it
        scattered from (source → flux), transposing the per-ℓ group-transfer
        :math:`\Sigma_{s,\ell}(g'\to g) \mapsto (g\to g')`.  Routes through
        :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_legendre_scattering_moments_transpose`
        — the SINGLE per-material dispatch site (Pattern 2), the transpose twin of
        :meth:`apply`'s ``apply_legendre_scattering_moments``.

        This is the bare Euclidean transpose; the metric-correct Hilbert adjoint
        :math:`\Lambda^{\dagger} = G^{-1}\Lambda^{T}G` is the
        :attr:`~orpheus.numerics.operator.LinearOperator.H` wrapper's job.  Because
        :math:`\Lambda` is the ONLY group-asymmetric factor of the frame-conjugated
        kernel :math:`R\circ\Lambda\circ M`, this method is what lets
        :math:`(R\circ\Lambda\circ M)^{T} = M^{T}\circ\Lambda^{T}\circ R^{T}` fall
        out of :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose`
        for free (campaign #276 A2).

        Typed :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
        in → :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
        out (the explicit role reversal); bare ndarray in → ndarray out (the
        endomorphic moment-space view the ``OperatorProduct.apply_transpose`` chain
        composes on).
        """
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux

        if isinstance(moments, HarmonicMomentSourceSink):
            out_values = self.mat_xs.apply_legendre_scattering_moments_transpose(
                moments.values, L=self.L, skip_l0=self.skip_l0,
            )
            return HarmonicMomentFlux.from_mesh_and_L(
                out_values, moments.mesh, moments.L,
                spatial_moments=moments.spatial_moments,
            )
        return self.mat_xs.apply_legendre_scattering_moments_transpose(
            moments, L=self.L, skip_l0=self.skip_l0,
        )

    @property
    def domain(self) -> "FunctionSpace":
        r"""The spherical-harmonic coefficient space — :math:`\Lambda` is endomorphic.

        :math:`\Lambda` acts block-diagonally per :math:`\ell` ON moment space, so its
        domain and codomain are BOTH the SH coefficient space
        :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace` of order :attr:`L`
        (``== frame.basis_space``). Carrying real spaces lets
        :class:`~orpheus.numerics.operator.OperatorProduct` admit the ``R∘Λ∘M``
        composition NATIVELY — the ``cast`` that papered over the former ``None``
        spaces (and short-circuited the composability guard) retires.
        """
        from orpheus.numerics.spaces import SphericalHarmonicSpace

        return SphericalHarmonicSpace.from_L(self.L)

    @property
    def codomain(self) -> "FunctionSpace":
        r"""Endomorphic on the SH coefficient space — codomain ``==`` :attr:`domain`."""
        return self.domain


@dataclass(frozen=True)
class N2NMomentOperator(LinearOperator):
    r"""The (n,2n) isotropic :math:`\ell=0` moment operator :math:`2\,\Sigma_{2n}`.

    The (n,2n) reaction is a DISTINCT isotropic (:math:`\ell=0`) group transfer —
    a *multiplication* channel (each event emits two neutrons), NOT scattering —
    so it is its own named operator rather than folded into
    :class:`LegendreMomentScattering`'s :math:`\ell=0` block. It is summed with
    :math:`\Lambda` in moment space (an
    :class:`~orpheus.numerics.operator.OperatorSum`) and the pair is
    frame-conjugated together, so the full isotropic + anisotropic in-scatter
    source is ONE :math:`R\circ(\Lambda + N_{2n})\circ M` (one analysis, one
    reconstruction). Campaign #276 A2 — physics-faithful: the multiplication
    reaction stays visible as a distinct operator with its own
    :meth:`apply` / :meth:`apply_transpose`, NOT hidden inside the scattering
    matmul.

    Endomorphic on the spherical-harmonic coefficient space of order :attr:`L`
    (it reads/writes ONLY the :math:`\ell=0` block, so it composes in an
    ``OperatorSum`` with :math:`\Lambda` on the same space); per-material dispatch
    lives in :meth:`MaterialXSField.apply_n2n_moments` (Pattern 2).
    """

    mat_xs: "MaterialXSField"
    L: int
    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )

    @property
    def is_adjointable(self) -> bool:
        # 2Σ_{2n}^T (apply_transpose) is the ℓ=0 group-transpose; caps ⊇
        # CAP_APPLY_TRANSPOSE. is_invertible inherits base False.
        return True

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r""":math:`2\,\Sigma_{2n}` applied to the :math:`\ell=0` moment (ℓ≥1 zero).

        Bare ndarray (the endomorphic moment-space view the ``frame.conjugate``
        ``OperatorProduct`` chain composes on — both for the forward
        :attr:`~ScatteringOperator.full_scatter_kernel` and its
        :meth:`apply_transpose`).
        """
        return self.mat_xs.apply_n2n_moments(moments)

    def apply_transpose(self, moments: np.ndarray) -> np.ndarray:
        r"""The :math:`\ell=0` group-transpose :math:`(2\,\Sigma_{2n})^{T}` (bare ndarray)."""
        return self.mat_xs.apply_n2n_moments_transpose(moments)

    @property
    def domain(self) -> "FunctionSpace":
        r"""The SH coefficient space of order :attr:`L` — :math:`N_{2n}` is endomorphic."""
        from orpheus.numerics.spaces import SphericalHarmonicSpace

        return SphericalHarmonicSpace.from_L(self.L)

    @property
    def codomain(self) -> "FunctionSpace":
        r"""Endomorphic — codomain ``==`` :attr:`domain`."""
        return self.domain


@dataclass
class ScatteringOperator(LinearOperator):
    r"""Scattering source operator :math:`S` (P0 + Pℓ + (n,2n)).

    Holds the precomputed per-material Legendre scattering matrices
    :math:`\Sigma_{s,\ell}`, the (n,2n) matrices :math:`\Sigma_{2n}`,
    the precomputed real spherical harmonics :math:`Y_\ell^m` evaluated
    at every quadrature direction, and a per-material cell-index map
    for vectorised assembly.

    Use :meth:`from_solver_data` to build instances from the same
    precomputed structures :class:`SNSolver` already holds (the canonical
    construction path used by Wave D Issue 13's bit-identical
    extraction).

    Attributes
    ----------
    mat_xs : MaterialXSField
        Macroscopic XS field — the single source of truth for both
        per-material scattering / (n,2n) data AND the cell-to-material
        topology (Issue #197 PR-TYPED-1).  Every per-material loop
        that previously lived on this class now routes through one of
        the typed verbs (``mat_xs.apply_p0_in_scatter``, ``apply_n2n``,
        etc.) which encapsulates the dispatch.
    quadrature : Quadrature
        The angular quadrature.  Carries ``N``, ``weights``, and
        :meth:`spherical_harmonics` — the previously-leaked
        ``n_ordinates`` / ``weights`` / ``Y`` constructor parameters
        all collapse into ``self.quadrature``.
    scattering_order : int
        Maximum Legendre order :math:`L` retained. ``0`` means P0 only.
    capabilities : frozenset[str]
        ``{"apply", "apply_transpose"}`` — :math:`S` has no efficient
        inverse (``solve``), but the adjoint :math:`S^{T}` is free via
        the harmonic-frame :attr:`full_scatter_kernel` (campaign #276
        A2b / #118): see :meth:`apply_transpose`.
    """

    mat_xs: "MaterialXSField"
    quadrature: "Quadrature"
    scattering_order: int

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )
    # Scattering is a BULK operator — the moment-folding `Σ_s · ⟨P_ℓ, ψ⟩`
    # reads and writes the bulk flux only (A_bb), no boundary action.
    # Issue #208 / Wave O. Class-level constant (unannotated so the
    # dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    # Lazy cache for the precomputed spherical harmonics — only
    # populated when ``scattering_order > 0`` (avoids paying the
    # cost on P0-only problems).
    _Y_cached: np.ndarray | None = field(
        default=None, init=False, repr=False,
    )

    #: The composite full-field space this operator acts on (P4.5 W-D).
    #: Threaded from the solver's ``sn_mesh.full_field_space`` via
    #: :meth:`from_solver_data`; ``None`` for the bare/test constructor (then
    #: ``domain``/``codomain`` report ``None`` and the composition guard skips
    #: — backward-compatible). When present it is the SAME ``(name, shape)``
    #: instance ``L``/``C``/``B`` carry, so the within-group ``(L+C) - S``
    #: :class:`~orpheus.numerics.operator.OperatorSum` guard validates the
    #: ``- S`` arm natively. ``S`` depends on this numerics ``FunctionSpace``,
    #: NOT on an SN mesh — the D5 cross-method-honest handle (relocation-ready
    #: for #261; ``S`` scatters in every method, not just SN).
    full_field_space: "FullFieldSpace | None" = field(
        default=None, repr=False, compare=False,
    )

    # ── Operator-algebra space metadata (P4.5 W-D) ───────────────────
    @property
    def domain(self) -> "FunctionSpace | None":
        r"""The composite full-field space (P4.5 W-D), or ``None`` if unthreaded.

        :math:`S` is a BULK 2-cell (the conjugation :math:`R\circ\Lambda\circ M`)
        but composes into the within-group loss ``(L+C) - S`` as a
        composite-field operator, so it advertises the SAME
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space` instance
        ``L``/``C``/``B`` carry (threaded via :meth:`from_solver_data`). Carrying
        the real space lets the residual
        :class:`~orpheus.numerics.operator.OperatorSum` guard VALIDATE the ``- S``
        arm instead of silently skipping it (the ``None``-spaced default
        pre-dates P4.5 W-D and keeps bare/test construction backward-compatible).
        :math:`S` reads and writes the bulk block only; the boundary block is
        inert, so domain == codomain on the composite.
        """
        return self.full_field_space

    @property
    def codomain(self) -> "FunctionSpace | None":
        # Endomorphic on the composite full-field space (see :meth:`domain`).
        return self.full_field_space

    @property
    def is_adjointable(self) -> bool:
        # S = R∘(Λ+N2N)∘M exposes its Euclidean transpose S^T via
        # :attr:`full_scatter_kernel` (the OperatorProduct adjoint chain);
        # caps ⊇ CAP_APPLY_TRANSPOSE. is_invertible inherits base False —
        # a scattering source operator is not invertible.
        return True

    # ── Convenience read-throughs (the legacy n_ordinates / nx / ny
    # / ng / weights / Y / cells_by_mat / sig_s / sig2 / sig_s0
    # attributes the test suite + internal call sites consumed).
    # Issue #197 PR-TYPED-1: these become read-throughs onto
    # mat_xs + quadrature so consumers don't break.  Marked TRANSIENT
    # in source comments where appropriate; full retirement deferred
    # to a follow-up PR once consumers all read mat_xs.* directly.

    @property
    def n_ordinates(self) -> int:
        """Number of angular ordinates :math:`N`."""
        return self.quadrature.N

    @property
    def spatial_shape(self) -> tuple[int, ...]:
        """Per-axis cell counts — read-through from :attr:`mat_xs`.

        C5.2 (#225): replaces the retired ``nx``/``ny`` pair (rank-
        generic; honest at any spatial dimension).
        """
        return self.mat_xs.spatial_shape

    @property
    def ng(self) -> int:
        """Energy group count."""
        return self.mat_xs.ng

    @property
    def weights(self) -> np.ndarray:
        """Quadrature weights ``(N,)``."""
        return self.quadrature.weights

    @property
    def Y(self) -> np.ndarray | None:
        r"""Real spherical harmonics :math:`Y_\ell^m(\Omega_n)`,
        shape ``(N, L+1, 2L+1)``, or ``None`` when ``L == 0``."""
        if self.scattering_order == 0:
            return None
        if self._Y_cached is None:
            self._Y_cached = self.quadrature.spherical_harmonics(
                self.scattering_order,
            )
        return self._Y_cached

    @cached_property
    def frame(self) -> HarmonicFrame:
        r"""The angular discrete frame — :class:`SphericalHarmonicBasis`
        (order :math:`L=` :attr:`scattering_order`) bound to the quadrature's
        :math:`S^2` measure.

        The SINGLE source of the analysis (:math:`M`) and reconstruction
        (:math:`R`) faces that realise the anisotropic Legendre
        redistribution :math:`R\circ\Lambda\circ M` (:attr:`kernel`). Both
        the §5.6 :attr:`kernel` and the angular-windowing in-sweep moment
        accumulation read THIS frame, so the projection table is shared
        term-for-term (no second evaluation, no risk of divergence —
        ``coding-elegance`` Pattern 2).

        A :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` — the
        carrier-typed specialization of the generic
        :class:`~orpheus.numerics.frame.GalerkinFrame`
        ``quadrature.angular_frame(L)`` (``from_galerkin`` reuses its basis +
        measure, so the inherited ndarray faces / :meth:`kernel` /
        :meth:`reconstruct_after` are bit-identical). The wrapper ADDS the typed
        verbs :meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.analyse`
        (:math:`M`) and
        :meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct`
        (:math:`R`) that the windowed in-scatter arm uses to reconstruct the
        per-ordinate source from the scattered moment SOURCE
        (:class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`)
        — the explicit, role-typed counterpart of the kernel's ndarray
        composition.

        ``frame.table`` is the :math:`Y_\ell^m(\hat\Omega_n)` tabulation,
        **bit-identical** to the legacy :attr:`Y` (both route through
        :meth:`SphericalHarmonicBasis.evaluate` on the same direction
        cosines). The :math:`S^2` embedding + binding lives once in
        :meth:`~orpheus.numerics.quadrature.Quadrature.angular_frame`.
        """
        return HarmonicFrame.from_galerkin(
            self.quadrature.angular_frame(self.scattering_order),
        )

    @property
    def sig_s(self) -> dict[int, list[np.ndarray]]:
        """TRANSIENT — per-material dense Legendre scattering dict.

        Routes through :meth:`MaterialXSField.sig_s_legendre`.  Kept
        as a shim for the four ``_build_rhs_*`` helpers in
        :mod:`orpheus.sn.solver` that still consume the per-material
        dict directly.  PR-TYPED-2 will rewire those helpers to
        consume ``mat_xs`` directly and retire this shim.
        """
        return {
            mid: self.mat_xs.sig_s_legendre(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig2(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material dense (n,2n) dict.  See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig_s0(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material P0 scattering matrix dict.
        See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.sig_s_legendre(mid)[0]
            for mid in self.mat_xs.materials
        }

    @property
    def cells_by_mat(self) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        """TRANSIENT — per-material cell-index dict.  See :attr:`sig_s`."""
        return self.mat_xs.cells_by_material

    # ── Anisotropic in-scatter: the moment→source map R·Λ_{ℓ≥1} ────────

    def _aniso_source_from_moment_values(
        self, moment_values: np.ndarray,
    ) -> np.ndarray:
        r"""Reconstruct the per-ordinate :math:`\ell\ge 1` in-scatter source
        from a flux-moment tensor — the **moment → source** half of the
        anisotropic scattering composition :math:`S_{\rm aniso}
        = \tfrac{1}{W}\,R\,\Lambda\,M`.

        .. math::

            (R\,\Lambda_{\ell\ge 1}\,\phi)_n(\vec r)

        where :math:`\Lambda_{\ell\ge 1}` is
        :class:`LegendreMomentScattering` (``skip_l0=True`` — the
        :math:`\ell=0` block is the P0 :meth:`add_iso_source` fast path)
        and :math:`R` is the :attr:`frame`'s
        :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face.
        The trailing :math:`1/W` is NOT applied here — that is the
        producer-side Pattern-7 normalisation at the :meth:`apply`
        boundary (lesson L18).

        It builds the :math:`R\,\Lambda` map as ``frame.reconstruct_after(Λ)`` —
        the §5.6 :attr:`kernel` sub-operator for inputs ALREADY in moment space
        (``coding-elegance`` Pattern 2: the composition IS the map). Its sole caller
        is the **windowed moment-iterate** :meth:`apply` arm (Phase 5a
        angular-windowing), whose iterate bulk IS the moment tensor :math:`\phi`, so
        :math:`M` is already done and only this :math:`R\,\Lambda` remains. The
        full-angular path (:meth:`build_aniso_source`) instead consumes the whole
        :math:`R\,\Lambda\,M` :attr:`kernel` (``= frame.conjugate(Λ)``); both share
        the frame's :math:`R` and :math:`\Lambda`, and the 0-ULP
        ``tests/sn/operators/test_scattering_kernel_crosscheck.py`` pins the
        ``kernel`` to that composition.

        It supersedes the per-:math:`\ell`
        ``_PerLegendreOrderScattering`` kernel (retired), which
        recomputed :math:`M\,\psi` independently for every Legendre
        order.  :math:`R` is linear, so reconstructing the full
        :math:`\Lambda`-scattered moment tensor at once equals summing the
        per-:math:`\ell` reconstructions (bit-identical in practice; pinned
        by the scattering regression tests).

        Parameters
        ----------
        moment_values : np.ndarray
            Flux-moment tensor ``(L+1, 2L+1, ng, nx, ny)`` — typically
            ``M.apply(psi)`` or the windowed iterate's
            :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
            ``.values``.

        Returns
        -------
        np.ndarray
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell\ge 1` in-scatter
            source, **pre** :math:`1/W`.
        """
        scatter = LegendreMomentScattering(
            mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True,
        )
        # R∘Λ as ONE typed operator (M already applied — the moments ARE M·ψ).
        return self.frame.reconstruct_after(scatter).apply(moment_values)

    @property
    def kernel(self) -> LinearOperator:
        r"""#257 S6 — the §5.6 integral kernel :math:`R \circ \Lambda \circ M`.

        Returns the anisotropic Legendre redistribution as a typed
        :class:`~orpheus.numerics.operator.OperatorProduct`

        .. math::

            \mathrm{kernel} \;=\; R \;\circ\; \Lambda_{\ell\ge 1}
                \;\circ\; M ,

        the genuinely-nonlocal-in-angle part of :math:`S`: it projects
        the per-ordinate flux onto harmonic moments (:math:`M`), applies
        the per-ℓ group-transfer cross sections on moment space
        (:math:`\Lambda_{\ell\ge 1}`, ``skip_l0=True``), and reconstructs
        the per-ordinate source via the addition theorem (:math:`R`). The
        factors are the frame faces and the moment-space scattering:

        * :math:`M` = ``frame.analysis`` (the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.analysis` face);
        * :math:`\Lambda_{\ell\ge 1}` = :class:`LegendreMomentScattering`
          ``(mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True)``;
        * :math:`R` = ``frame.reconstruction`` (the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face);

        composed via ``frame.conjugate(Λ)`` into
        ``OperatorProduct(R, OperatorProduct(Λ, M))`` (whose
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply` is
        ``a.apply(b.apply(x))``), giving
        :math:`\mathrm{kernel.apply}(\psi.\mathrm{values}) =
        R(\Lambda(M\,\psi))`. The production aniso paths CONSUME this operator:
        :meth:`build_aniso_source` is ``(1/W)·kernel``; the windowed
        moment-iterate arm consumes the sub-operator
        ``frame.reconstruct_after(Λ)`` (:math:`R\,\Lambda`).

        With this property :class:`ScatteringOperator` satisfies the
        §5.6
        :class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
        Protocol: scattering IS a nonlocal integral kernel (the ``R∘Λ∘M``
        angular redistribution), while the isotropic
        :math:`\ell=0` :math:`P_0` in-scatter fast path
        (:meth:`add_iso_source`) and the :math:`(n,2n)` doubling
        (:meth:`add_n2n_source`) are the LOCAL / separate components of
        the full :meth:`apply`.

        The production path (Cardinal Rule 2 — the composition IS the realization)
        --------------------------------------------------------------------------

        ``kernel`` is the production :math:`\ell\ge 1` map, NOT a parallel
        "semantic" reading of a hand-chain. :meth:`build_aniso_source` IS
        ``(1/W)·kernel`` (it calls ``self.kernel.apply(ψ)``, then the producer-side
        :math:`1/W`), and the windowed moment-iterate :meth:`apply` arm calls the
        sub-operator ``frame.reconstruct_after(Λ)`` (:math:`R\,\Lambda` — :math:`M`
        already done). There is ONE composition now; the 0-ULP
        ``tests/sn/operators/test_scattering_kernel_crosscheck.py`` is its
        DEFINITIONAL identity (it was a twin-guard between two parallel chains).
        The producer-side :math:`1/W` normalisation (lesson L18) lives OUTSIDE this
        kernel, at the :meth:`apply` boundary, so :math:`\mathrm{kernel.apply}`
        returns the per-ordinate :math:`\ell\ge 1` source **pre**-:math:`1/W`.

        Built fresh on each access (the factors are cheap reference bindings) to
        honor the :attr:`mat_xs` / :attr:`Y` read-through.

        Raises
        ------
        ValueError
            If ``scattering_order == 0`` — an isotropic operator has no
            anisotropic integral kernel (the :math:`\ell\ge 1` redistribution is
            empty; :attr:`Y` is ``None``). The §5.6 Kernel surface is meaningful
            only for anisotropic scattering; the P0 in-scatter is the LOCAL
            component handled by :meth:`add_iso_source`.

        Returns
        -------
        LinearOperator
            ``frame.conjugate(Λ)`` ``= OperatorProduct(R, OperatorProduct(Λ, M))``
            — the typed ``R∘Λ∘M`` anisotropic redistribution.
        """
        if self.scattering_order == 0:
            raise ValueError(
                "ScatteringOperator.kernel requires scattering_order >= 1: "
                "an isotropic (P0-only) operator has no anisotropic integral "
                "kernel (the R∘Λ∘M angular redistribution is empty). The P0 "
                "in-scatter is the LOCAL component, handled by add_iso_source."
            )
        # R ∘ Λ ∘ M as ONE typed operator: the frame conjugates Λ (per-ℓ
        # moment-space scattering) between its analysis (M) and reconstruction (R)
        # faces. Λ now carries real spaces (== frame.basis_space), so the
        # OperatorProduct composability guard validates R∘Λ∘M natively — NO cast.
        return self.frame.conjugate(
            LegendreMomentScattering(
                mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True,
            )
        )

    @property
    def full_scatter_kernel(self) -> LinearOperator:
        r"""The FULL in-scatter source kernel :math:`R\circ(\Lambda_{\ell\ge 0} + N_{2n})\circ M`.

        The modernized production scattering source (campaign #276 A2a): the
        COMPLETE P0 + anisotropic + (n,2n) in-scatter as ONE frame-conjugated
        operator. The isotropic ℓ=0 scattering (:math:`\Lambda_{\ell=0}`), the
        anisotropic ℓ≥1 redistribution (:math:`\Lambda_{\ell\ge 1}`), and the
        (n,2n) multiplication (the DISTINCT :class:`N2NMomentOperator`) are summed
        in moment space (an :class:`~orpheus.numerics.operator.OperatorSum`) and
        conjugated by the angular frame TOGETHER — one analysis :math:`M`, one
        reconstruction :math:`R`. The per-ordinate scattering source is
        ``(1/W)·full_scatter_kernel.apply(ψ)`` (this is what :meth:`apply`
        computes); its transpose is ``(1/W)·full_scatter_kernel.apply_transpose(ψ*)``
        (the adjoint scattering :math:`S^{T}`, campaign #276 A2b).

        This is the iso path "made nice like aniso": the isotropic in-scatter now
        rides the SAME frame conjugation as the anisotropic redistribution, so the
        whole operator's transpose falls out of
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose` for free
        (the M/R face transposes + :math:`\Lambda^{T}` + :math:`N_{2n}^{T}`).

        Distinct from :attr:`kernel` (the §5.6 ℓ≥1 ANISOTROPIC subcomponent — the
        :class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
        protocol surface + the 0-ULP crosscheck oracle): this is the FULL
        ℓ≥0 + (n,2n) source. Built fresh per access (read-through to :attr:`mat_xs`,
        so depletion / thermal-feedback updates show up immediately).
        """
        return self.frame.conjugate(
            LegendreMomentScattering(
                mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=False,
            )
            + N2NMomentOperator(mat_xs=self.mat_xs, L=self.scattering_order)
        )

    @cached_property
    def isotropic_kernel(self) -> LinearOperator:
        r"""The model-shared isotropic in-scatter energy operator
        :math:`\Sigma_{s,0} + 2\,\Sigma_{2n}` on the scalar flux.

        The :math:`\ell=0` energy source — P0 in-scatter plus the :math:`(n,2n)`
        doubling — as the cross-model
        :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
        ``+``
        :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
        :class:`~orpheus.numerics.operator.OperatorSum`.  This is the SAME
        per-material ``mat_xs`` dispatch the inline :meth:`add_iso_source` /
        :meth:`add_n2n_source` verbs used, so
        :meth:`_assemble_per_ordinate_source` routes the SN forward isotropic
        source through it **bit-identically** (campaign #276 P2): the
        :class:`~orpheus.numerics.operator.OperatorSum`'s
        ``a.apply(φ) + b.apply(φ)`` equals the in-place ``add_iso`` then
        ``add_n2n`` accumulation (pinned 0-ULP by
        ``test_isotropic_scattering.py::test_sum_equals_inplace_iso_then_n2n``).

        The same energy operation lives — re-implemented — in every transport
        model (CP / MoC / diffusion / homogeneous / MC); these operators are the
        model-independent extraction (``coding-elegance`` Pattern 2 over the six
        solver families).  Space-anonymous (``space=None``): the bulk
        scalar-flux energy action carries no composition-guard space at this
        producer-side use.  Cached: the kernel reads the XS *through*
        :attr:`mat_xs` at ``apply`` time, so depletion / thermal-feedback updates
        to the same field show up immediately (the operator snapshots no values).
        """
        from orpheus.transport.operators.isotropic_scattering import (
            IsotropicN2N,
            IsotropicScattering,
        )
        return IsotropicScattering(self.mat_xs) + IsotropicN2N(self.mat_xs)

    def _assemble_per_ordinate_source(
        self,
        phi: "np.ndarray | ScalarFlux",
        aniso_or_none: "AngularSourceSink | None",
        mesh: "SNMesh",
    ) -> "AngularSourceSink":
        r"""Combine the P0 + (n,2n) iso source (from the scalar flux
        :math:`\phi_0`) with the pre-:math:`/W` :math:`\ell\ge 1` aniso source
        into the per-ordinate scattering source :math:`(\text{iso}/W) +
        \text{aniso}`.

        The single source of truth (``coding-elegance`` Pattern 2) for the
        producer-side :math:`1/W` combine (Pattern 7 / lesson L18) shared by
        every ``apply`` arm that emits a per-ordinate
        :class:`~orpheus.transport.source_sinks.AngularSourceSink` — the
        full-angular :class:`AngularFlux` arm and the windowed
        :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
        arm.  Both reduce to "compute :math:`\phi` and the pre-:math:`/W`
        aniso, then assemble"; the :math:`/W` convention lives HERE, once.

        Parameters
        ----------
        phi : np.ndarray or ScalarFlux
            Scalar flux :math:`\phi_0` (iso magnitude) driving P0 + (n,2n).
        aniso_or_none : AngularSourceSink or None
            Per-ordinate :math:`\ell\ge 1` source ALREADY in per-ordinate
            magnitude (post-:math:`/W`), or ``None`` for
            ``scattering_order == 0``.
        mesh : SNMesh
            Phase-space carrier (sizes the zero accumulators + ``sum_w``).
        """
        # The iso source carries the trailing spatial-moment axis when the
        # driving flux does (the φ̂ iterate, #240 D5b-S3 — the Σ_s ⊗ I lift
        # scatters every spatial moment).  Read the width OFF the ScalarFlux's
        # space (the single source of truth); a bare-ndarray φ stays scalar.
        spatial_moments = _spatial_moments_of(phi)
        # Route the isotropic (P0 + (n,2n)) energy source through the
        # model-shared K_iso operators (#276 P2): IsotropicScattering = Σ_s0,
        # IsotropicN2N = 2Σ_2n, summed.  SAME per-material ``mat_xs`` dispatch as
        # the legacy add_iso/add_n2n verbs ⇒ bit-identical (OperatorSum.apply ≡
        # the in-place accumulation, pinned 0-ULP by test_isotropic_scattering).
        iso: ScalarSourceSink = ScalarSourceSink.from_mesh(
            self.isotropic_kernel.apply(phi), mesh, spatial_moments=spatial_moments,
        )
        aniso = (
            aniso_or_none
            if aniso_or_none is not None
            else AngularSourceSink.zeros_on(mesh, spatial_moments=spatial_moments)
        )
        sum_w = float(self.quadrature.weights.sum())
        return (iso / sum_w) + aniso

    @classmethod
    def from_solver_data(
        cls,
        *,
        mat_xs: "MaterialXSField",
        quadrature: "Quadrature",
        scattering_order: int,
        full_field_space: "FullFieldSpace | None" = None,
    ) -> "ScatteringOperator":
        """Construct from a :class:`MaterialXSField` + quadrature.

        ``full_field_space`` (P4.5 W-D) is the composite
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space` the solver
        threads so ``S``'s ``domain``/``codomain`` match ``L``/``C``/``B``
        and the within-group ``(L+C) - S`` composition guard validates the
        ``- S`` arm. ``None`` (the default) leaves ``S`` space-less — the
        guard skips it, preserving the legacy contract for direct callers.

        Issue #197 PR-TYPED-1 — the constructor surface collapses the
        eight separate per-material handles (``sig_s``, ``sig2``,
        ``sig_s0``, ``cells_by_mat``, ``Y``, ``weights``, ``n_ordinates``,
        ``nx``, ``ny``, ``ng``) into two typed objects.  The
        :class:`MaterialXSField` carries everything per-material plus
        the spatial topology; the :class:`Quadrature` carries
        ``N`` / ``weights`` / harmonics.
        """
        return cls(
            mat_xs=mat_xs,
            quadrature=quadrature,
            scattering_order=scattering_order,
            full_field_space=full_field_space,
        )

    # ── In-place helpers (preserve bit-identity vs SNSolver pre-Wave-D) ─

    def add_iso_source(
        self,
        Q: "np.ndarray | ScalarSourceSink",
        phi: "np.ndarray | ScalarFlux",
    ) -> "np.ndarray | ScalarSourceSink | None":
        r"""Add P0 in-scatter source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += sig_s0[mid].T @ phi[:, ic, jc]`` where
        ``sig_s0[mid]`` is ``(ng, ng)`` indexed ``[g_from, g_to]``.

        Issue #197 PR-TYPED-3 introduces the typed-action overload:

        * **Raw-in (legacy)** — ``Q: np.ndarray`` is mutated in place
          and ``None`` is returned (the Wave A–D / PR-TYPED-1 contract).
        * **Typed-in (return-new)** — ``Q: ScalarSourceSink`` is treated
          as an immutable input; a NEW :class:`ScalarSourceSink` is
          returned carrying ``Q.values + in_scatter`` (Pattern 4 —
          frozen typed inputs stay immutable; Pattern 1 — the caller
          spells the algebra as ``Q = scattering.add_iso_source(Q,
          phi)``).

        Parameters
        ----------
        Q : np.ndarray or ScalarSourceSink
            Isotropic source carrier.  Raw ``(ng, nx, ny)`` ndarray is
            mutated in place; typed :class:`ScalarSourceSink` returns a
            new instance.
        phi : np.ndarray or ScalarFlux
            Scalar flux.  Either form is accepted; the underlying
            values are unwrapped before the per-material dispatch.

        Returns
        -------
        np.ndarray or ScalarSourceSink or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`ScalarSourceSink`.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material dispatch lives ONLY
        inside :meth:`MaterialXSField.apply_p0_in_scatter`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.source_sinks import ScalarSourceSink
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSourceSink):
            Q_values = Q.values.copy()
            self.mat_xs.apply_p0_in_scatter(Q_values, phi_values)
            # Preserve Q's spatial-moment width (#240 D5b-S3 — the φ̂ accumulator
            # carries the trailing 2^d axis; re-wrapping without it would lose
            # the typed factor and raise on the shape).
            return ScalarSourceSink.from_mesh(
                Q_values, Q.mesh, spatial_moments=Q.spatial_moments_per_axis,
            )
        self.mat_xs.apply_p0_in_scatter(Q, phi_values)
        return None

    def add_n2n_source(
        self,
        Q: "np.ndarray | ScalarSourceSink",
        phi: "np.ndarray | ScalarFlux",
    ) -> "np.ndarray | ScalarSourceSink | None":
        r"""Add (n,2n) source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += 2 · sig2[mid].T @ phi[:, ic, jc]``.

        Issue #197 PR-TYPED-3 introduces the same typed-action overload
        as :meth:`add_iso_source` — raw-in mutates in place and returns
        ``None``; typed-in returns a fresh :class:`ScalarSourceSink`.

        Parameters
        ----------
        Q : np.ndarray or ScalarSourceSink
            Isotropic source carrier.
        phi : np.ndarray or ScalarFlux
            Scalar flux.

        Returns
        -------
        np.ndarray or ScalarSourceSink or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`ScalarSourceSink`.

        Notes
        -----
        Issue #197 PR-TYPED-1 — per-material dispatch lives ONLY
        inside :meth:`MaterialXSField.apply_n2n`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.source_sinks import ScalarSourceSink
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSourceSink):
            Q_values = Q.values.copy()
            self.mat_xs.apply_n2n(Q_values, phi_values)
            # Preserve Q's spatial-moment width (#240 D5b-S3, as add_iso_source).
            return ScalarSourceSink.from_mesh(
                Q_values, Q.mesh, spatial_moments=Q.spatial_moments_per_axis,
            )
        self.mat_xs.apply_n2n(Q, phi_values)
        return None

    def build_aniso_source(
        self,
        angular_flux: "np.ndarray | AngularFlux | None",
    ) -> "np.ndarray | AngularSourceSink | None":
        r"""Build per-ordinate Pℓ (:math:`\ell \ge 1`) scattering source.

        Implements the Galerkin reconstruction :eq:`pn-scatter` from
        the angular-flux moments :eq:`flux-moments` (declared in
        :ref:`theory-discrete-ordinates`) as the **literal operator
        composition**

        .. math::

            Q^{\rm aniso}_n(\vec r) \;=\;
            \frac{1}{W}\,(R \, \Lambda \, M \, \psi)_n(\vec r)

        of the three primitives shipped in
        :mod:`orpheus.numerics.projection` and
        :class:`LegendreMomentScattering` (Wave 1 / C1.1):

        * :math:`M` — the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.analysis` face:
          :math:`\psi(N, \cdot) \mapsto \phi^{\ell m}(L+1, 2L+1, \cdot)`
          via :math:`\phi_\ell^m = \sum_n w_n Y_\ell^m(\Omega_n) \psi_n`.
        * :math:`\Lambda` :class:`LegendreMomentScattering`: per-ℓ
          block-diagonal cross-section action on moment space (skip ℓ=0,
          handled by P0 in-scatter).
        * :math:`R` — the :attr:`frame`'s
          :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face:
          :math:`\phi^{\ell m} \mapsto q_n` via the addition-theorem
          reconstruction :math:`q_n = \sum_\ell (2\ell+1) \sum_m
          Y_\ell^m(\Omega_n) \phi_\ell^m`.

        The :math:`(2\ell+1)` factor in :math:`R` is the addition-theorem
        scaling sourced from
        :attr:`~orpheus.numerics.spaces.SphericalHarmonicSpace.addition_theorem_factor`
        — the single canonical home of the SH convention's
        :math:`(2\ell+1)` literal (ERR-039 endpoint).  :math:`R` is the
        addition-theorem reconstruction, NOT the W-weighted Hilbert
        adjoint of :math:`M` (the two differ by the
        :class:`SphericalHarmonicSpace` Gram :math:`g_C`); the
        relationship is pinned by
        ``tests/numerics/test_spherical_harmonic_space.py``.

        The trailing :math:`1/W = 1/\sum_n w_n` is the **producer-side
        per-ordinate projection** introduced in R-1 Step 4 A1: the Pℓ
        aniso source enters the per-ordinate transport equation
        :math:`(\Omega\cdot\nabla + \sigma_t)\psi_n = Q/W +
        q_{\rm aniso,n}` already in per-ordinate magnitude, so the
        sweep does NOT need to apply ``/W`` again
        (cf. :func:`~orpheus.sn.loss_representation.transport_sweep` whose
        ``aniso_source`` parameter is documented as per-ordinate
        magnitude post-A1).  See ``coding-elegance`` SKILL.md
        §"Convention crosswalk template" / lesson L18.

        Total flop count is identical to the legacy hand-rolled
        ``for n in range(N)`` triple loop; the iteration over ordinates
        is now internal to :func:`numpy.einsum` inside :math:`R` and
        :math:`M`, **not** a Python loop.

        Parameters
        ----------
        angular_flux : np.ndarray or AngularFlux or None
            Angular flux shape ``(N, ng, nx, ny)`` (Issue #196
            PR-INDEX-4 — principled) from the most recent sweep, or
            ``None`` on the first source iteration before any sweep
            has been run.  Issue #197 PR-TYPED-3 — typed
            :class:`AngularFlux` is accepted; the result is then a
            typed :class:`AngularSourceSink` (preserving the type chain
            through the scattering composition).

        Returns
        -------
        np.ndarray or AngularSourceSink or None
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell \ge 1`
            contribution in **per-ordinate magnitude** (the trailing
            ``/W`` is applied here, R-1 Step 4 A1).  Returns ``None``
            when ``scattering_order == 0`` or ``angular_flux is None``.
            Type matches the input.

        Notes
        -----
        The :attr:`frame`'s analysis and reconstruction faces are
        layout-agnostic in their trailing axes — only the leading ordinate /
        harmonic axes are consumed by the einsums.  Switching ψ from ``(N, nx, ny,
        ng)`` to ``(N, ng, nx, ny)`` therefore transparently produces a
        moment field shape ``(L+1, 2L+1, ng, nx, ny)`` (instead of
        ``(L+1, 2L+1, nx, ny, ng)``), which is exactly what the
        principled :class:`LegendreMomentScattering` consumes after
        PR-INDEX-4.
        """
        if self.scattering_order == 0 or angular_flux is None:
            return None
        # D-I.2 (2026-05-29): only the typed :class:`AngularFlux` is
        # accepted.  The bare-ndarray fallthrough retired alongside the
        # singledispatch ``np.ndarray`` arm on :meth:`apply`; the typed
        # leaf carries the mesh, which is required to construct the
        # output :class:`AngularSourceSink`.
        mesh = angular_flux.mesh
        # S_aniso = (1/W)·kernel.  The §5.6 :attr:`kernel` (= ``frame.conjugate(Λ)``,
        # the typed R∘Λ∘M redistribution) is THE production composition; this arm
        # adds only the producer-side :math:`/sum_w` (Pattern 7 / lesson L18, which
        # lives OUTSIDE the kernel).  The windowed moment-iterate arm consumes the
        # sub-operator ``frame.reconstruct_after(Λ)`` (M already done) — same R, Λ.
        sum_w = float(self.weights.sum())
        out_values = self.kernel.apply(angular_flux.values) / sum_w
        # The R·Λ·M chain is spatial-moment-axis-agnostic (#240 D5b-S3); when the
        # angular iterate carries φ̂, the aniso source does too — the typed wrap
        # selects the SpatialMomentSpace factor (read off the iterate's space).
        return AngularSourceSink.from_mesh(
            out_values, mesh,
            spatial_moments=angular_flux.spatial_moments_per_axis,
        )

    # ── Foldable / residual split ─────────────────────────────────────
    #
    # The Phase G four-operator unification splits scattering as
    #
    #     S = S_foldable + S_residual
    #
    # where ``S_foldable`` is the **P0 within-group self-scatter**
    # only: per ordinate the source is :math:`\Sigma_{s,0}^{g\to g}
    # \phi_g`, direction-independent, so the cell-balance denominator
    # can absorb it as :math:`\sigma_r = \sigma_t -
    # \Sigma_{s,0}^{g\to g}` (the "removal" cross-section). Every
    # other piece of :math:`S` lives in ``S_residual``:
    #
    # * **Cross-group P0** (off-diagonal of ``sig_s[mid][0]``) — couples
    #   distinct energy groups, cannot collapse into a per-cell scalar.
    # * **All :math:`P_\ell \ge 1`** — direction-dependent through
    #   :math:`Y_\ell^m(\Omega_n)`, fundamentally unfoldable into
    #   :math:`\sigma_r`. This is NOT an SN-specific quirk; it is
    #   inherent to anisotropic scattering's angular structure for
    #   every transport solver family (SN, CP, MoC, diffusion-P_N).
    # * **(n,2n) doubling** — emits two neutrons per absorption;
    #   folding into a "removal" cross-section :math:`\sigma_r` is
    #   conceptually wrong even when ``sig2[mid][g, g]`` is non-zero.
    #
    # Both halves are :class:`ScatteringOperator` siblings (same class,
    # different per-material arrays). No new wrapper type. The split
    # is purely additive: ``S.apply(\psi) =
    # S.foldable_part().apply(\psi) + S.residual_part().apply(\psi)``
    # at FP-non-associativity precision (typically ``rtol=1e-14``).
    #
    # The split is the correct DATA decomposition of the within-group
    # self-scatter — and it is exactly the input a Diffusion Synthetic
    # Acceleration (DSA, issue #2) preconditioner consumes (Σ_s0 → the
    # diffusion removal coefficient). This substep lands the data API only;
    # no solver/sweep/iteration code consumes the new methods yet.
    #
    # ⚠ LATENT CORRECTNESS TRAP (issue #215 — measured 2026-06-05): do NOT
    # wire the foreseen ``ψ_{n+1} = A_wg.solve(S_residual·ψ + q)`` with a
    # σ_r-SWEEP as ``A_wg.inverse()``.  The σ_r-sweep inverts
    # ``(Ω·∇ + σ_r·I)`` — a DIAGONAL-in-angle removal — whereas
    # ``S_foldable = Σ_s0·P_iso`` is the ISOTROPIC-PROJECTION self-scatter
    # operator (rank-1 in angle, ``φ/Σw``).  The two coincide ONLY when ψ is
    # isotropic.  Wiring ``A_wg.solve(S_residual)`` with the σ_r-sweep is
    # therefore exact on a fully-reflective uniform box (isotropic flux) but
    # ships **46–56 % flux errors on any anisotropic problem** (vacuum /
    # heterogeneous) — it passes the isotropic unit tests and silently
    # corrupts real cases.  The exact variant (keep the ``−Σ_s0⊙ψ`` remnant
    # on the RHS) has the right fixed point but DIVERGES (``Σ_s0/σ_r ≈ 39``).
    # The stable + correct within-group self-scatter fold is CONSISTENT DSA
    # (#2) or KRYLOV (already production, splitting-invariant, rate-optimal
    # on every BC).  Any future within-group accelerator MUST be gated on an
    # ANISOTROPIC config — the isotropic box cannot see this error.

    def foldable_part(self) -> "ScatteringOperator":
        r"""Return the P0 within-group self-scatter sibling of :math:`S`.

        Carries only the diagonal of ``sig_s[mid][0]`` per material —
        the within-group self-scatter cross-section
        :math:`\Sigma_{s,0}^{g\to g}` per energy group. All other
        scattering channels (cross-group P0, all :math:`P_\ell \ge 1`,
        and (n,2n)) live in :meth:`residual_part`.

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order = 0``, diagonal-only
            ``sig_s[mid][0]``, zero ``sig2``, and a derived
            :class:`MaterialXSField` carrying the foldable overrides.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material loop that built the
        ``sig_s_foldable`` / ``sig2_foldable`` dicts lives ONLY inside
        :meth:`MaterialXSField.foldable_sig_s`.
        """
        sig_s_foldable = self.mat_xs.foldable_sig_s()
        sig2_foldable = {
            mid: np.zeros_like(self.mat_xs.n2n_matrix(mid))
            for mid in self.mat_xs.materials
        }
        derived_mat_xs = self.mat_xs.with_overridden_sig_s_and_n2n(
            sig_s_dense=sig_s_foldable,
            n2n_dense=sig2_foldable,
        )
        return ScatteringOperator(
            mat_xs=derived_mat_xs,
            quadrature=self.quadrature,
            scattering_order=0,
            full_field_space=self.full_field_space,
        )

    def residual_part(self) -> "ScatteringOperator":
        r"""Return the non-foldable sibling of :math:`S`.

        Carries everything :meth:`foldable_part` does not: the
        off-diagonal of ``sig_s[mid][0]`` (cross-group P0), every
        :math:`P_\ell \ge 1` block verbatim, and ``sig2[mid]``
        verbatim.

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order ==
            self.scattering_order``, zero-diagonal P0 matrix, and
            unchanged :math:`P_\ell \ge 1` + (n,2n) data.

        Notes
        -----
        Algebraic contract:
        ``S.apply(\psi) \approx S.foldable_part().apply(\psi) +
        S.residual_part().apply(\psi)`` at ``rtol=1e-14``.

        Issue #197 PR-TYPED-1 — the per-material loop that built the
        ``sig_s_residual`` dict lives ONLY inside
        :meth:`MaterialXSField.residual_sig_s`.
        """
        sig_s_residual = self.mat_xs.residual_sig_s()
        # (n,2n) carried verbatim — pull from mat_xs.
        sig2_residual = {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }
        derived_mat_xs = self.mat_xs.with_overridden_sig_s_and_n2n(
            sig_s_dense=sig_s_residual,
            n2n_dense=sig2_residual,
        )
        return ScatteringOperator(
            mat_xs=derived_mat_xs,
            quadrature=self.quadrature,
            scattering_order=self.scattering_order,
            full_field_space=self.full_field_space,
        )

    def is_foldable_into_sigma_r(self) -> bool:
        r"""Return ``True`` iff this operator is structurally the
        :meth:`foldable_part` of some parent :math:`S`.

        Returns
        -------
        bool
            ``True`` iff ``scattering_order == 0`` AND every material's
            P0 is diagonal AND every material's (n,2n) is zero.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material allclose checks live
        ONLY inside :meth:`MaterialXSField.is_p0_diagonal_with_zero_n2n`.
        """
        if self.scattering_order != 0:
            return False
        return self.mat_xs.is_p0_diagonal_with_zero_n2n()

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""Return the per-material foldable cross-section :math:`(\sigma_{s,0}^{g\to g})_g`.

        Returns
        -------
        dict[int, np.ndarray]
            ``{mid: (ng,) array}`` keyed by material id. Each array
            is a fresh copy.

        Notes
        -----
        Issue #197 PR-TYPED-1 — delegates to
        :meth:`MaterialXSField.foldable_sigma`.
        """
        return self.mat_xs.foldable_sigma()

    # ── LinearOperator surface ─────────────────────────────────────────

    @singledispatchmethod
    def _apply_impl(self, psi) -> "Any":
        r"""Runtime dispatch for :meth:`apply` — see the typed overloads.

        Applies the full scattering source operator :math:`S\,\psi`.
        Dispatched on the input *carrier* type via
        :func:`functools.singledispatchmethod`.  The public :meth:`apply`
        name aliases this dispatcher at runtime and carries the honest
        per-carrier ``@overload`` typing surface (#257 S8c), so callers
        statically see the exact output type for each input carrier:

        * :class:`~orpheus.transport.timed_full_field.TimedFullField`
          → :class:`~orpheus.transport.full_field.FullField` — composite
          bulk + boundary variant.  Bulk follows the full :math:`P_\ell`
          Galerkin path; boundary is the implicit-zero
          :class:`~orpheus.transport.source_sinks.BoundarySourceSink`
          (Option β3 — scattering is volumetric; Wave O #208 will encode
          the bulk-only nature in the type).  #257 S8a made the codomain
          the timeless :class:`FullField` (the matvec leaf is a base
          arrow; the iteration driver reattaches the timed type).
        * :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
          → :class:`~orpheus.transport.source_sinks.ScalarSourceSink` —
          :math:`P_0` in-scatter + :math:`(n,2n)` doubling only, in **iso
          scalar magnitude**.  No :math:`P_\ell` (scalar flux has no
          angular info); no :math:`1/W` (scalar consumers — diffusion, CP,
          kinetics — do not project to per-ordinate).
        * :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
          → :class:`~orpheus.transport.source_sinks.AngularSourceSink` —
          full :math:`P_\ell` Galerkin reconstruction in **per-ordinate
          magnitude** (the trailing :math:`1/W` projection lives at this
          producer boundary; Pattern 7).
        * :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
          → :class:`~orpheus.transport.source_sinks.AngularSourceSink` —
          the Phase 5a angular-windowing path: :math:`S` consumes flux
          MOMENTS (the moments ARE :math:`M\psi`, so the :math:`M`
          projection is skipped), bit-identical to the
          :class:`AngularFlux` arm for :math:`\phi = M\psi`.

        The single :class:`ScatteringOperator` instance serves every
        consumer kind via type-directed dispatch — Pattern 1 (read-as-the-
        math) + Pattern 7 (producer-side normalisation, the :math:`1/W`
        lives at the apply boundary, decided by consumer type).  See
        ``coding-elegance`` SKILL.md §"Convention crosswalk template" and
        lesson L18.

        Notes
        -----
        Internal helpers :meth:`add_iso_source`, :meth:`add_n2n_source`,
        and :meth:`build_aniso_source` remain available for callers
        that need the iso / aniso pieces separately.
        :meth:`build_aniso_source` returns **per-ordinate magnitude**
        (the :math:`1/W` is applied there); the per-ordinate
        :class:`AngularFlux` variant of :meth:`apply` combines the iso
        piece (in iso magnitude, projected by :math:`1/W` at the
        boundary) with the aniso piece (already per-ordinate) into a
        single :class:`AngularSourceSink`.
        """
        raise TypeError(
            f"ScatteringOperator.apply: unsupported input type "
            f"{type(psi).__name__}; expected TimedFullField, ScalarFlux, "
            f"AngularFlux, or HarmonicMomentFlux.  Dispatch table is "
            f"registered via @singledispatchmethod."
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        r"""Composite :class:`FullField` variant — bulk-only scattering source.

        Registered on the timeless :class:`FullField` (W-C): a
        :class:`TimedFullField` iterate dispatches here via MRO (it IS a
        ``FullField``), so the runtime is behaviour-preserving, and a bare
        ``FullField`` now dispatches correctly. Reads only ``psi.bulk``
        (history-blind). The ``@overload`` static stubs name ``FullField``
        too (W-F), matching this runtime registration.

        Math: identical to the :class:`AngularFlux` branch above —
        reduce ``psi.bulk`` angular → scalar, build iso :math:`P_0 +
        (n,2n)` source on the typed
        :class:`~orpheus.transport.source_sinks.ScalarSourceSink` accumulator,
        build the per-ordinate :math:`P_\ell\ge 1` Galerkin
        contribution via :meth:`build_aniso_source` (now accepting
        the L2 :class:`AngularFlux` on the composite's bulk), then
        combine via the producer-side :math:`1/W` projection
        (Pattern 7).  The output bulk is a pure-Field
        :class:`AngularFlux`; the output boundary is the
        **implicit-zero** :class:`BoundaryFlux` — scattering is
        volumetric, no face-trace contribution.

        Per Option β3 (`#208
        <https://github.com/deOliveira-R/ORPHEUS/issues/208>`_) the
        implicit-zero boundary is a transitional shim: Wave O will
        introduce :class:`BulkOperator` / :class:`FullOperator`
        Protocols so that scattering's volumetric nature is encoded
        in the *type*, not in a zero-valued boundary member.  Until
        then the composite return enables ``L.apply(state) -
        S.apply(state) - F.apply(state)`` to compose under
        :meth:`TimedFullField.__sub__` once all four operators expose
        the composite branch (D-H.1c).
        """
        # Delegate the bulk source to the bulk-type dispatch arm and wrap
        # with the implicit-zero boundary.  ``psi.bulk`` is either the
        # full-angular :class:`AngularFlux` (1-D / curvilinear / un-windowed
        # 2-D) OR — Phase 5a angular-windowing — the reduced-moment
        # :class:`HarmonicMomentFlux` (the 2-D Cartesian windowed SI
        # iterate); both bulk arms return the same per-ordinate
        # :class:`AngularSourceSink`, so this composite arm is one
        # delegation (Pattern 2 — no duplicated iso/n2n/aniso assembly).
        # B.5.2: the operator output IS a source (Sψ rate density).  The
        # boundary is the implicit-zero :class:`BoundarySourceSink` —
        # scattering is volumetric (Option β3 / Wave O #208).
        # #257 S8c: ``psi.bulk`` is statically the broad ``BulkField``; the
        # cast names the runtime truth documented above (AngularFlux OR
        # HarmonicMomentFlux — both dispatch arms return AngularSourceSink)
        # so the typed ``apply`` overloads resolve.
        combined = self.apply(cast("AngularFlux | HarmonicMomentFlux", psi.bulk))
        # #257 S8a: history-free (the matvec leaf is a base arrow
        # ``FullField -> FullField``; the comonad lives on the driver, which
        # reattaches the timed type when this source is added to the timed rhs).
        return FullField(
            bulk=combined,
            boundary=BoundarySourceSink.zeros_on(cast("SNMesh", psi.bulk.mesh)),
        )

    @_apply_impl.register
    def _(self, phi: ScalarFlux) -> "ScalarSourceSink":
        r"""Typed ScalarFlux variant — iso scalar magnitude output (P0 + n2n only).

        Math:
        :math:`Q(\vec r, g) = \Sigma_{s,0}(g'\to g)\,\phi_{g'}(\vec r)
        + 2\,\Sigma_{2n}(g'\to g)\,\phi_{g'}(\vec r)`.
        No :math:`P_\ell` (scalar flux lacks the angular info needed
        for the Galerkin reconstruction); no :math:`1/W` (consumers in
        scalar-flux equations — diffusion, CP, kinetics outer — do not
        project to per-ordinate).

        **Deliberately retained (W-F, 2026-06-26) — a named-future-consumer
        surface, NOT dead weight.**  This arm has no current production
        caller: the within-group SI/Krylov path feeds the composite
        :class:`FullField` / per-ordinate :class:`AngularFlux` /
        :class:`HarmonicMomentFlux` arms, never a bare :class:`ScalarFlux`,
        and (unlike fission's internally-reached :class:`ScalarFlux` arm) no
        sibling arm delegates here — it is orphan at both ends today.  It is
        kept as the typed entry-point for the scalar-carrier cross-method
        consumers (diffusion / CP / kinetics outer) that the cross-method
        field architecture (`#205
        <https://github.com/deOliveira-R/ORPHEUS/issues/205>`_) will wire.
        Per the coding-standards retirement audit this keep-vs-retire call is
        recorded EXPLICITLY rather than left a silent orphan (user steered
        keep, 2026-06-26).
        """
        mesh = phi.mesh
        iso: ScalarSourceSink = ScalarSourceSink.zeros_on(mesh)
        iso = self.add_iso_source(iso, phi)
        iso = self.add_n2n_source(iso, phi)
        return iso

    @_apply_impl.register
    def _(self, psi: AngularFlux) -> "AngularSourceSink":
        r"""Typed :class:`AngularFlux` variant — per-ordinate magnitude output.

        Math: identical to the :class:`TimedFullField` arm above on the
        bulk axis — reduce ``psi`` angular → scalar, build iso :math:`P_0
        + (n,2n)` source on the typed
        :class:`~orpheus.transport.source_sinks.ScalarSourceSink` accumulator,
        build the per-ordinate :math:`P_\ell\ge 1` Galerkin contribution
        via :meth:`build_aniso_source`, then combine via the producer-side
        :math:`1/W` projection (Pattern 7).

        D-I.2 (2026-05-29) added this typed leaf alongside the existing
        :class:`TimedFullField` composite arm; the legacy bare-ndarray
        :class:`numpy.ndarray` arm retired in the same commit.  Direct
        :class:`AngularFlux` consumers (the within-group iteration when
        scattering is computed on the bulk-only carrier) now have a
        first-class dispatch arm instead of routing through a transient
        :class:`TimedFullField` wrap.
        """
        # φ = ∫ψ dΩ (scalar), aniso = (1/W) R Λ M ψ (per-ordinate), then the
        # shared producer-side assembly. The iso (P0 + n2n) keeps the cheap
        # reaction-rate fast path (NO moment tensor) — it is a load-bearing
        # PERF optimisation on the SI-sweep hot path; routing it through the
        # frame regresses LD/P0 badly (campaign #276 A2a finding). The frame
        # form lives on as :attr:`full_scatter_kernel` for the (non-hot-path)
        # adjoint transpose only (A2b).
        return self._assemble_per_ordinate_source(
            psi.integrate_angular(), self.build_aniso_source(psi), psi.mesh,
        )

    @_apply_impl.register
    def _(self, phi_moments: HarmonicMomentFlux) -> "AngularSourceSink":
        r"""Windowed moment-iterate variant — :math:`S` consumes flux MOMENTS.

        The Phase 5a angular-windowing path.  When the within-group SI
        iterate is stored as harmonic moments :math:`\phi_\ell^m` (the 2-D
        Cartesian windowed iterate) instead of the full per-ordinate
        :class:`AngularFlux`, :math:`S` consumes the moments WITHOUT the
        :math:`M` projection — the moments ARE :math:`M\psi`, so projecting
        again would be redundant work.

        Structurally parallel to the :class:`AngularFlux` arm, and
        **bit-identical** to it for :math:`\phi = M\psi` (de-risk proven):

        * the :math:`\ell=0` block IS the scalar flux — ORPHEUS uses
          unnormalized real harmonics (:math:`Y_0^0 = 1`), so
          :meth:`HarmonicMomentFlux.scalar_flux` equals
          :meth:`AngularFlux.integrate_angular` term-for-term (the typed
          accessor single-sources that convention); it feeds the identical
          P0 + (n,2n) fast path;
        * the :math:`\ell\ge 1` aniso takes the EXPLICIT typed grid path:
          :math:`\Lambda` scatters the flux moments to the in-scatter source
          moments (``HarmonicMomentFlux → HarmonicMomentSourceSink`` — the
          role-changing edge), then the frame's :math:`R`
          (:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.reconstruct`)
          reconstructs the per-ordinate :class:`AngularSourceSink`. This equals
          the full-angular path's :math:`R\,\Lambda` after its own :math:`M`,
          and the kernel's ndarray ``reconstruct_after(Λ)`` reference
          (:meth:`_aniso_source_from_moment_values`, the 0-ULP crosscheck's
          oracle) — same :math:`\Lambda` kernel + :math:`R` face — but with the
          scattered moment source materialized as a TYPED carrier, so the
          flux → source → angular role flow reads off the signatures.

        Output convention matches the :class:`AngularFlux` arm exactly: both
        end at the shared :meth:`_assemble_per_ordinate_source` (per-ordinate
        :class:`AngularSourceSink`, producer-side :math:`1/W`, Pattern 7).
        """
        mesh = phi_moments.mesh
        if self.scattering_order == 0:
            aniso = None
        else:
            # The explicit typed grid path: Λ scatters the flux moments to the
            # in-scatter SOURCE moments (HarmonicMomentFlux → HarmonicMomentSourceSink
            # — the role-changing edge), then the frame's R reconstructs the
            # per-ordinate AngularSourceSink (HarmonicMomentSourceSink → AngularSourceSink
            # — role-preserving), then the producer-side 1/W (Pattern 7). The
            # scattered moment source is materialized as a TYPED carrier, so the
            # flux→source→angular role flow reads off the signatures rather than
            # an opaque ndarray. Numerically equals the kernel's ndarray
            # reconstruct_after(Λ) reference (same Λ kernel + R face); the
            # spatial-moment φ̂ axis threads through the carriers automatically.
            scattered = LegendreMomentScattering(
                mat_xs=self.mat_xs, L=self.scattering_order, skip_l0=True,
            ).apply(phi_moments)
            aniso = self.frame.reconstruct(scattered) / float(self.weights.sum())
        # ℓ=0 moment IS the scalar flux (Y_0^0 = 1) — the typed accessor
        # carries that convention (== integrate_angular bit-exactly).
        return self._assemble_per_ordinate_source(
            phi_moments.scalar_flux(), aniso, mesh,
        )

    if TYPE_CHECKING:
        # Honest per-carrier typing surface (#257 S8c).  ``ScatteringOperator``
        # is NOT an endomorphism ``V -> V`` (the mixin's nominal contract):
        # it maps each input carrier to a DISTINCT output carrier.  These
        # ``@overload`` stubs exist only for the type checker; the public
        # ``apply`` IS the runtime dispatcher (``apply = _apply_impl`` below),
        # so callers statically see e.g. ``S.apply(ScalarFlux) ->
        # ScalarSourceSink`` instead of the dispatcher's untyped fallback.
        @overload
        def apply(self, psi: FullField, /) -> "FullField": ...
        @overload
        def apply(self, phi: ScalarFlux, /) -> "ScalarSourceSink": ...
        @overload
        def apply(self, psi: AngularFlux, /) -> "AngularSourceSink": ...
        @overload
        def apply(
            self, phi_moments: HarmonicMomentFlux, /,
        ) -> "AngularSourceSink": ...
        def apply(self, x: Any, /) -> Any: ...
    else:
        apply = _apply_impl

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r"""The adjoint scattering source :math:`S^{T}\chi =
        (1/W)\,\mathrm{full\_scatter\_kernel}^{T}\chi` (campaign #276 A2b, closes
        `#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_).

        :math:`S` maps a flux to the per-ordinate in-scatter source it emits;
        :math:`S^{T}` is the group-and-angle transpose the adjoint transport
        equation :math:`(L+C-S)^{T}\psi^{*}=q^{*}` needs.  Because the production
        FORWARD keeps the scalar fast-path (:attr:`isotropic_kernel` +
        :meth:`build_aniso_source`) for SI-sweep performance (the A2a regression),
        the adjoint — NOT the hot path — rides the validated harmonic-frame form
        instead: :attr:`full_scatter_kernel` :math:`= R\circ(\Lambda_{\ell\ge0}
        + N_{2n})\circ M`, whose transpose :math:`M^{T}\circ(\Lambda+N_{2n})^{T}
        \circ R^{T}` falls out of
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply_transpose` for
        free (the Phase-D :math:`M/R` face transposes :math:`+\,\Lambda^{T}+
        N_{2n}^{T}`).  ONE expression gives the COMPLETE transpose — iso
        :math:`\ell=0` + aniso :math:`\ell\ge1` + (n,2n); the producer-side
        :math:`1/W` (Pattern 7) transposes as the scalar it is
        (:math:`(A/W)^{T}=A^{T}/W`).

        This is the **Euclidean** transpose (the plain group-and-angle matvec
        adjoint, L12) — NOT the metric Hilbert adjoint ``.H`` (which would carry
        the angular Gram).  It is pinned by the reciprocity
        :math:`\langle S\psi,\chi\rangle=\langle\psi,S^{T}\chi\rangle` against the
        INDEPENDENT forward fast-path
        (``test_scattering_adjoint.py::TestFullScatterKernel``): a genuine
        cross-check, since :meth:`apply` (scalar fast-path) and
        :meth:`apply_transpose` (frame form) are structurally different
        representations of the same operator.  The forward/transpose asymmetry is
        principled and oracle-pinned — the forward equivalence
        :math:`(1/W)\,\mathrm{full\_scatter\_kernel}.\mathrm{apply}\equiv
        S.\mathrm{apply}` holds to ~1e-12
        (``test_reproduces_forward_scattering_source``).

        Parameters
        ----------
        chi : np.ndarray or carrier
            The daggered per-ordinate field :math:`\chi` of shape
            ``(N, ng, *spatial)``; a typed carrier's ``.values`` are unwrapped.
            Returns a bare ``ndarray`` — the typed adjoint-source carrier
            (``AdjointFlux`` / ``Solution.adjoint``) arrives in #116 (A5).
        """
        chi_values = np.asarray(getattr(chi, "values", chi))
        return self.full_scatter_kernel.apply_transpose(chi_values) / float(
            self.weights.sum()
        )

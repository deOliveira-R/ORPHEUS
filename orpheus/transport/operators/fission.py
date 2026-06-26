r"""Multigroup fission source operator as a :class:`LinearOperator`.

This module owns the **fission source operator** :math:`F` from the
operator-algebra view of the Boltzmann transport equation,

.. math::

    (L - S - F)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (L - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`L` is the streaming + collision operator, :math:`S` is
the scattering operator (see :mod:`orpheus.transport.operators.scattering`), and
:math:`F` is the **fission emission operator**

.. math::

    (F\,\phi)_g(\vec r) = \chi_g\,
    \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r).

The structure is **rank-1 in energy**: the production rate
:math:`\sum_{g'} \nu\Sigma_{f,g'}\phi_{g'}` is a scalar per cell and
the emission spectrum :math:`\chi_g` redistributes it across groups via
an outer product. This rank-1 factorisation is the algebraic origin of
why power iteration converges geometrically in a critical reactor
(dominance ratio :math:`|k_1/k_0|`) — the eigenvalue problem
:math:`A^{-1}F\phi = k\phi` has :math:`F` of rank
:math:`O(N_{\rm cells})` per group, but only one **per cell** in
energy because :math:`\chi` is a rank-1 spectrum.

The §5.6 Kernel reading — fission as an integral kernel
=======================================================

In the grand-report §5.6 suffix law (see
:mod:`orpheus.transport.operators.integral_kernel_operator`) fission is a
**Kernel**: a *nonlocal* operator whose action integrates the flux
against a measure on the group axis (the emission at :math:`(\vec r, g)`
reads the flux at *every* group :math:`g'`). #257 S6 makes this a
type-level fact — :class:`FissionOperator` exposes both the integral
structure :attr:`~FissionOperator.kernel` (the rank-1 ``χ ⊗ νΣf``
:class:`~orpheus.numerics.operator.TensorProductOperator`, present since
Wave T) and the §5.6 middle factor :attr:`~FissionOperator.production_rate`
(the S5 :class:`~orpheus.transport.production_rate_functional.ProductionRateFunctional`),
and therefore satisfies the
:class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
Protocol.

The **semantic** decomposition (Frame 3 of the cross-domain-attacker
memo) is the composition

.. math::

    F \;=\; M_\chi \;\circ\; \mathrm{ProductionRate} \;\circ\;
        M_{\nu\Sigma_f} ,

where :math:`M_{\nu\Sigma_f}` multiplies the flux by the production
cross section, :math:`\mathrm{ProductionRate}` (the
:attr:`~FissionOperator.production_rate` functional) contracts the group
axis to the per-cell density :math:`p(\vec r)`, and :math:`M_\chi`
re-broadcasts that density across groups weighted by :math:`\chi`. This
is the *reading*, NOT the realization: the fused
:class:`~orpheus.numerics.operator.RankOneOperator`-backed
:attr:`~FissionOperator.kernel` is the efficient evaluation
(``coding-elegance`` Pattern 5 — the right primitive, not the unfolded
product). The two agree byte-for-byte
(:math:`\chi \cdot \mathrm{production\_rate.evaluate}(\phi) \equiv
F.\mathrm{apply}(\phi)` at 0 ULP) because
:meth:`~orpheus.transport.production_rate_functional.ProductionRateFunctional.evaluate`
reproduces the rank-1 ``inner`` reduction — a fact **pinned** by
``tests/sn/operators/test_fission_kernel_crosscheck.py`` (B.2), not
merely asserted here. The matvec arms are UNCHANGED in S6 (additive,
bit-identical).

Per Cardinal Rule 2 (architecture) this lifts the
``SNSolver.compute_fission_source`` math out of the solver and into a
single operator object. The math is **moved verbatim** (Wave D Issue 13
is a bit-identical extraction). The eigenvalue division by :math:`k`
**stays at the solver level** — :meth:`apply` returns
:math:`F\,\phi`, NOT :math:`F\,\phi / k`. The caller (the EigenvalueSolver
Protocol's ``compute_fission_source``, which is :class:`SNSolver`'s
delegator) divides by :math:`k`. Two reasons:

1. The :class:`LinearOperator` Protocol contract is *linear*:
   :meth:`apply` returns :math:`F\,\phi`, not :math:`F\,\phi/k`. The
   :math:`1/k` factor is a scalar multiple in the algebra
   :math:`F^{eff} = (1/k)\,F` — express it via :class:`ScaledOperator`
   if the algebra needs it explicitly.

2. The eigenvalue iteration owns :math:`k`. The fission operator
   should not be re-built every outer iteration just because :math:`k`
   changed; the operator state is :math:`(\chi, \nu\Sigma_f)`, which
   is fixed across the outer iteration.

Capability advertisement
========================

:pydata:`capabilities = frozenset({CAP_APPLY})`. No ``solve``: the
operator is **rank-1 in energy** per cell — its inverse does not exist.
No ``apply_transpose``: the adjoint fission operator is
:math:`F^T\,\phi^* = \nu\Sigma_f \cdot (\chi \cdot \phi^*)`, structurally
distinct (it sums over the **adjoint** energy distribution); the
current ORPHEUS solver does not consume :math:`F^T`. Both can be added
in a future wave when the adjoint transport machinery lands.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import singledispatchmethod
from typing import TYPE_CHECKING, Any, overload

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    CAP_APPLY,
    IdentityOperator,
    LinearOperator,
    RankOneOperator,
    TensorProductOperator,
)

# Runtime imports for :func:`singledispatchmethod.register` — see
# ``scattering.py`` for the same pattern.  These types form a leaf in
# the SN dependency graph (they do not import fission.py).  The L2
# pure-Field :class:`AngularFlux` is re-aliased to ``AngularFlux`` to
# disambiguate from the legacy ``orpheus.sn.angular_flux.AngularFlux``
# which still rides on the operator-algebra path until D-H.1c.
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.production_rate_functional import ProductionRateFunctional
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import ScalarSourceSink
from orpheus.transport.timed_full_field import TimedFullField

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces import FullFieldSpace


__all__ = ["FissionOperator"]


@dataclass
class FissionOperator(LinearOperator):
    r"""Fission source operator :math:`F = \chi\,\otimes\,\nu\Sigma_f`.

    Reads :math:`\chi(\vec r)` and :math:`\nu\Sigma_{f,g}(\vec r)`
    through a :class:`MaterialXSField` (Issue #197 PR-TYPED-1) — the
    same per-cell typed views every other operator (L, C, S)
    consumes.  The action is a contraction over groups (the production
    rate) followed by a broadcast across the emission spectrum.

    Use :meth:`from_solver_data` to build instances; pass
    ``mat_xs=sn_mesh.material_xs_field()``.

    Attributes
    ----------
    mat_xs : MaterialXSField
        Macroscopic XS field carrying ``emission_spectrum`` (χ) and
        ``fission_production`` (νΣ_f) per-cell views.
    capabilities : frozenset[str]
        ``{"apply"}`` — the rank-1 structure forbids a useful inverse,
        and the adjoint surface is not yet a consumer.
    """

    mat_xs: "MaterialXSField"

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    #: The composite full-field space (P4.5 W-D) — threaded from the solver's
    #: ``sn_mesh.full_field_space`` via :meth:`from_solver_data`. ``F`` NEVER
    #: enters a production :class:`~orpheus.numerics.operator.OperatorSum`
    #: (it is always applied separately as ``F.apply(ψ)`` — the fission source
    #: ``Fψ/k``), so this is the honest cross-method tag, NOT a guard-activator:
    #: it keeps ``F`` uniform with ``C``/``S`` and ready for a future
    #: ``(L+C-S-F-B)`` posing / DSA composition. ``None`` ⟹ ``domain``/
    #: ``codomain`` report ``None`` (guard skips; backward-compatible). ``F``
    #: depends on this numerics ``FunctionSpace``, NOT an SN mesh (D5).
    full_field_space: "FullFieldSpace | None" = field(
        default=None, repr=False, compare=False,
    )

    # Fission is a BULK operator — `χ · νΣ_f · ⟨1, ψ⟩` reads and writes
    # the bulk flux only (A_bb), no boundary action. Issue #208 / Wave O.
    # Class-level constant (unannotated so the dataclass does not treat
    # it as a field).
    block_role = BlockRole.BULK

    # ── Read-through properties for backwards compatibility ────────────

    @property
    def chi(self) -> np.ndarray:
        r""":math:`\chi(\vec r)` per-cell view, shape ``(ng, nx, ny)``.

        Read-through onto :attr:`mat_xs.emission_spectrum`.
        """
        return self.mat_xs.emission_spectrum

    @property
    def sig_p(self) -> np.ndarray:
        r""":math:`\nu\Sigma_f(\vec r)` per-cell view, shape ``(ng, nx, ny)``.

        Read-through onto :attr:`mat_xs.fission_production`.
        """
        return self.mat_xs.fission_production

    # ── Operator-algebra space metadata (P4.5 W-D) ───────────────────
    @property
    def domain(self) -> "FunctionSpace | None":
        r"""The composite full-field space (P4.5 W-D), or ``None`` if unthreaded.

        :math:`F` advertises the SAME
        :attr:`~orpheus.sn.geometry.SNMesh.full_field_space` instance
        ``L``/``C``/``S``/``B`` carry (threaded via :meth:`from_solver_data`).
        Unlike :math:`C`/:math:`S`, :math:`F` never enters a production
        :class:`~orpheus.numerics.operator.OperatorSum` (the fission source is
        applied as ``F.apply(ψ)`` and divided by ``k`` at the eigenvalue
        layer), so this space activates no production guard today — it is the
        honest cross-method tag that makes a future ``(L+C-S-F-B)`` posing or
        DSA composition type-checkable. :math:`F` reads the bulk block only;
        domain == codomain on the composite.
        """
        return self.full_field_space

    @property
    def codomain(self) -> "FunctionSpace | None":
        # Endomorphic on the composite full-field space (see :meth:`domain`).
        return self.full_field_space

    @classmethod
    def from_solver_data(
        cls, *, mat_xs: "MaterialXSField",
        full_field_space: "FullFieldSpace | None" = None,
    ) -> "FissionOperator":
        """Construct from a :class:`MaterialXSField`.

        Issue #197 PR-TYPED-1 — the constructor surface collapses the
        ``(chi, sig_p)`` ndarray pair into one :class:`MaterialXSField`
        handle that carries both views consistently with the rest of
        the four-operator algebra. ``full_field_space`` (P4.5 W-D) is the
        composite space the solver threads so ``F.domain``/``codomain`` match
        the rest of the algebra (the honest cross-method tag; ``F`` never
        composes in production, so ``None`` is harmless for direct callers).
        """
        return cls(mat_xs=mat_xs, full_field_space=full_field_space)

    @property
    def kernel(self) -> "TensorProductOperator":
        r"""Wave T step T.2 — the rank-1 TP kernel of the fission action.

        Returns the 2-factor :class:`TensorProductOperator`

        .. math::

            \mathrm{RankOneOperator}(\chi,\,\nu\Sigma_f,\,\mathrm{axis}=0)
            \;\&\; \mathrm{IdentityOperator}()

        which composes :math:`F = \chi \otimes \nu\Sigma_f` (Grand
        Report v3 §15 line 2008) as a type-visible separable form
        per §16A.10's ``B = G_patch \otimes K_omega \otimes K_g``
        decomposition.  The first factor encodes the rank-1
        group-axis outer product (with spatial parameters carried
        in the ``(ng, nx, ny)``-shaped ``left`` and ``right``); the
        second factor advertises the spatial-axis broadcast.

        Built fresh on each access to honor the read-through semantics
        for :attr:`mat_xs`: depletion or thermal-feedback updates that
        mutate :attr:`MaterialXSField.emission_spectrum` or
        :attr:`MaterialXSField.fission_production` in-place show up
        immediately in the next :meth:`apply` call.  The construction
        is O(1) — two array-reference bindings plus a constructor
        call — so the lookup cost is negligible compared to the
        ``np.einsum`` inside :meth:`RankOneOperator.apply`.

        Bit-identity with the legacy two-step formulation
        ``np.einsum("gxy,gxy->xy", sig_p, phi) * chi[None, :, :]``
        is preserved because :meth:`RankOneOperator.apply` performs
        the same ``(right * x).sum(axis=0, keepdims=True)`` reduction
        followed by ``left * inner`` broadcast — the IEEE-754
        pairwise reduction order is identical.

        Returns
        -------
        TensorProductOperator
            The 2-factor TP whose first element is a
            :class:`RankOneOperator(self.chi, self.sig_p, axis=0)`
            and whose second element is an :class:`IdentityOperator`.

        Notes
        -----
        The verification gate for Wave T step T.2 inspects this
        property — the type signature change from the legacy
        opaque-einsum body to a typed :class:`TensorProductOperator`
        is the deliverable.  See
        ``tests/sn/test_fission_operator.py::TestRankOneTensorProductKernel``.
        """
        return (
            RankOneOperator(self.chi, self.sig_p, axis=0)
            & IdentityOperator()
        )

    @property
    def production_rate(self) -> ProductionRateFunctional:
        r"""#257 S6 — the §5.6 middle factor :math:`\mathrm{ProductionRate}`.

        Returns the S5
        :class:`~orpheus.transport.production_rate_functional.ProductionRateFunctional`
        over this operator's production cross section
        :math:`\nu\Sigma_f` (read through the S2 typed accessor
        :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.fission_production_field`).
        It contracts the group axis of a flux to the per-cell fission
        emission density

        .. math::

            p(\vec r) \;=\; \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,
                \phi_{g'}(\vec r) ,

        the **middle** factor of the §5.6 semantic decomposition
        :math:`F = M_\chi \circ \mathrm{ProductionRate} \circ
        M_{\nu\Sigma_f}` (Frame 3 of the cross-domain-attacker memo
        ``coefficient_field_promotion_frames.md``). Naming the
        contraction turns the most physically central diagnostic in
        criticality — the per-cell fission rate — from an anonymous
        einsum into a typed, inspectable
        :class:`~orpheus.numerics.functional.Functional`
        (``coding-elegance`` Pattern 3).

        Together with :attr:`kernel` (present since Wave T) this property
        makes :class:`FissionOperator` satisfy the §5.6
        :class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
        Protocol — it carries the operator surface (``apply`` /
        ``capabilities``) AND the integral-kernel structure.

        Bit-identity (the realization is the fused
        :class:`~orpheus.numerics.operator.RankOneOperator`)
        --------------------------------------------------------

        ``production_rate`` NAMES the §5.6 *semantic* middle factor; it
        does NOT replace the matvec. The realization stays the fused
        :attr:`kernel` (a
        :class:`~orpheus.numerics.operator.RankOneOperator`-backed
        :class:`~orpheus.numerics.operator.TensorProductOperator`), per
        ``coding-elegance`` Pattern 5 (the right primitive, not the
        unfolded product). The two AGREE byte-for-byte by construction:
        :meth:`ProductionRateFunctional.evaluate` reproduces the
        ``inner = (νΣf · φ).sum(axis=0, keepdims=True)`` reduction inside
        :meth:`~orpheus.numerics.operator.RankOneOperator.apply` (same
        numpy primitive, same axis, same ``keepdims``), and the χ
        re-broadcast reproduces ``RankOneOperator``'s ``left * inner`` —
        so :math:`\chi \cdot \mathrm{production\_rate.evaluate}(\phi)
        \equiv F.\mathrm{apply}(\phi)` at 0 ULP. This is pinned by
        ``tests/sn/operators/test_fission_kernel_crosscheck.py`` (B.2).

        Built fresh on each access (O(1) — one
        :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
        wrap) to honor the :attr:`mat_xs` read-through: a depletion /
        thermal-feedback update to :math:`\nu\Sigma_f` shows up
        immediately in the next access.
        """
        return ProductionRateFunctional(self.mat_xs.fission_production_field)

    @singledispatchmethod
    def _apply_impl(self, phi) -> "Any":
        r"""Runtime dispatch for :meth:`apply` — see the typed overloads.

        Applies :math:`F\,\phi`: emission rate × spectrum, no :math:`1/k`.

        .. math::

            (F\,\phi)_g(\vec r) = \chi_g(\vec r)\,
              \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

        Dispatched on the input *carrier* type via
        :func:`functools.singledispatchmethod`.  The public :meth:`apply`
        name aliases this dispatcher at runtime and carries the honest
        per-carrier ``@overload`` typing surface (#257 S8c), so callers
        statically see the exact output type for each input carrier:

        * :class:`~orpheus.transport.timed_full_field.TimedFullField`
          → :class:`~orpheus.transport.full_field.FullField` — composite
          bulk + boundary variant.  The bulk is the iso fission source
          projected to per-ordinate magnitude; the boundary is the
          implicit-zero :class:`~orpheus.transport.source_sinks.BoundarySourceSink`
          (Option β3 — fission has no boundary action; Wave O #208 will
          encode the bulk-only nature in the type).  #257 S8a made the
          codomain the timeless :class:`FullField` (the matvec leaf is a
          base arrow; the iteration driver reattaches the timed type).
        * :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
          → :class:`~orpheus.transport.source_sinks.ScalarSourceSink` —
          fission emission in **iso scalar magnitude**, for consumers in
          scalar-flux equations (eigenvalue outer / depletion / diffusion)
          that do not project to per-ordinate.  Symmetric with
          :meth:`ScatteringOperator.apply`'s :class:`ScalarFlux` variant.
        * :class:`numpy.ndarray` — bare ``(ng, *spatial)`` iso scalar
          fission source.  Preserved for ``KEigenvalue`` / depletion /
          diffusion consumers that feed bare arrays at outer-iteration
          boundaries.

        Fission has no :math:`P_\ell` aniso component (the emission
        spectrum is isotropic by construction); see ``coding-elegance``
        SKILL.md §"Convention crosswalk template" and lesson L18 for the
        Pattern 7 producer-side normalisation discipline.
        """
        raise TypeError(
            f"FissionOperator.apply: unsupported input type "
            f"{type(phi).__name__}; expected TimedFullField, ScalarFlux, "
            f"or numpy.ndarray.  Dispatch table is registered via "
            f"@singledispatchmethod."
        )

    @_apply_impl.register
    def _(self, psi: FullField) -> "FullField":
        r"""Composite :class:`FullField` variant — bulk-only fission emission.

        Registered on the timeless :class:`FullField` (W-C): a
        :class:`TimedFullField` iterate dispatches here via MRO (it IS a
        ``FullField``), so the runtime is behaviour-preserving, and a bare
        ``FullField`` now dispatches correctly. Reads only ``psi.bulk``
        (history-blind). The ``@overload`` static stubs still name
        ``TimedFullField`` — restructuring them is W-F (the overload retire).

        Math: identical to the :class:`AngularFlux` branch above —
        reduce bulk angular → scalar via :math:`\phi = \sum_n w_n
        \psi_n`, compute iso fission source :math:`F\phi`, project to
        per-ordinate via :math:`/W`.  The output bulk is a pure-Field
        :class:`AngularFlux`; the output boundary is an **implicit
        zero** :class:`BoundaryFlux` because the fission operator
        has no boundary action (the emission spectrum
        :math:`\chi(\vec r)` lives only on cell-centred volumes).

        Per Option β3 (`#208
        <https://github.com/deOliveira-R/ORPHEUS/issues/208>`_) the
        implicit-zero boundary is a transitional shim: Wave O will
        introduce :class:`BulkOperator` /
        :class:`FullOperator` Protocols so that fission's bulk-only
        nature is encoded in the *type*, not in a zero-valued boundary
        member.  Until then the composite return enables
        :math:`L.\mathrm{apply}(\psi) - S.\mathrm{apply}(\psi) -
        F.\mathrm{apply}(\psi)` to compose under
        :meth:`TimedFullField.__sub__` once all four operators expose
        the composite branch (D-H.1c).
        """
        phi_scalar = psi.bulk.integrate_angular()
        # Reuse the ScalarFlux branch — single source of truth for the
        # per-cell production-rate × emission-spectrum contraction.
        fission_iso: ScalarSourceSink = self.apply(phi_scalar)
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            BoundarySourceSink,
        )
        per_ord = AngularSourceSink.from_isotropic(
            fission_iso.values, psi.bulk.mesh,
        )
        return FullField(
            # B.5.2: the operator output IS a source (Fψ rate density) — emit
            # the AngularSourceSink directly, not a re-wrap into AngularFlux.
            # #257 S8a: history-free (the matvec leaf is a base arrow; the
            # comonad lives on the driver, which reattaches the timed type when
            # this source is added to the timed rhs).
            bulk=per_ord,
            boundary=BoundarySourceSink.zeros_on(psi.bulk.mesh),
        )

    @_apply_impl.register
    def _(self, phi: ScalarFlux) -> "ScalarSourceSink":
        r"""Typed ScalarFlux variant — iso scalar magnitude output.

        Math:
        :math:`Q_g(\vec r) = \chi_g(\vec r)\,
        \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)`.
        Iso scalar magnitude — no :math:`1/W` (scalar consumers do
        not project).

        Returns :class:`ScalarSourceSink` (R-1 Step 4 A1 — was
        :class:`ScalarFlux` pre-A1; the return type now matches
        :meth:`ScatteringOperator.apply`'s ScalarFlux variant by
        symmetry, and reflects the dimensional truth that the
        fission output is a source quantity, not a flux).

        Wave T step T.2: the inner math is delegated to the typed
        :attr:`kernel` (a 2-factor :class:`TensorProductOperator`
        wrapping :class:`RankOneOperator`).  This dispatch arm
        handles only the typed-flux layer (extraction of
        ``phi.values`` and packaging into :class:`ScalarSourceSink`);
        the rank-1 outer-product math itself lives at the L1
        primitive level.
        """
        out = self.kernel.apply(phi.values)
        return ScalarSourceSink.from_mesh(out, phi.mesh)

    @_apply_impl.register
    def _(self, phi_arr: np.ndarray) -> np.ndarray:
        r"""Bare-ndarray legacy variant — iso scalar fission source.

        Shape contract: input ``(ng, nx, ny)`` scalar flux, output
        ``(ng, nx, ny)`` iso fission source.  Preserved for
        ``KEigenvalue`` / depletion / diffusion outer-iteration
        consumers that still feed bare ``(ng, nx, ny)`` arrays.  No
        type wrapping; the bare path bypasses the type layer entirely.

        Wave T step T.2: delegates to :attr:`kernel` (same as the
        typed :class:`ScalarFlux` branch).  Single source of truth
        for the rank-1 outer-product math.
        """
        return self.kernel.apply(phi_arr)

    if TYPE_CHECKING:
        # Honest per-carrier typing surface (#257 S8c).  ``FissionOperator``
        # is NOT an endomorphism ``V -> V`` (the mixin's nominal contract):
        # it maps each input carrier to a DISTINCT output carrier.  These
        # ``@overload`` stubs exist only for the type checker; the public
        # ``apply`` IS the runtime dispatcher (``apply = _apply_impl``
        # below), so callers statically see e.g. ``F.apply(ScalarFlux) ->
        # ScalarSourceSink`` instead of the dispatcher's untyped fallback.
        @overload
        def apply(self, phi: TimedFullField, /) -> "FullField": ...
        @overload
        def apply(self, phi: ScalarFlux, /) -> "ScalarSourceSink": ...
        @overload
        def apply(self, phi: np.ndarray, /) -> np.ndarray: ...
        def apply(self, phi: Any, /) -> Any: ...
    else:
        apply = _apply_impl

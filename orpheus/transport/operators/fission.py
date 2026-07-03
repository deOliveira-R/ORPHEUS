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

The §5.6 Kernel reading — fission as a rank-1 integral kernel
============================================================

In the grand-report §5.6 suffix law (see
:mod:`orpheus.transport.operators.integral_kernel_operator`) fission is a
**Kernel**: a *nonlocal* operator whose action integrates the flux
against a measure on the group axis (the emission at :math:`(\vec r, g)`
reads the flux at *every* group :math:`g'`). :class:`FissionOperator`
exposes that integral structure as :attr:`~FissionOperator.kernel` and
therefore satisfies the
:class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
Protocol.

Fission is the **rank-1 dyad** :math:`F = |\chi\rangle\langle\nu\Sigma_f|`:

.. math::

    F \;=\; \texttt{outer}\bigl(\chi,\
        \mathrm{ReactionRateFunctional}(\nu\Sigma_f)\bigr) ,

a reconstruction column :math:`\chi` (the emission spectrum) tensored
with the production-rate row co-vector :math:`\langle\nu\Sigma_f| =
\mathrm{ReactionRateFunctional}(\nu\Sigma_f)`. The matvec **routes
through** the functional: :meth:`apply` is :math:`\chi \cdot
\mathrm{production\_rate.evaluate}(\phi) = \chi \cdot \langle\nu\Sigma_f,
\phi\rangle`, so the production-rate co-vector IS the contraction the
kernel performs — there is no separate "fused" realization to drift from
the named factor (the procedural twin the earlier S5 design carried is
**dissolved**). This is the rank-1 (single-mode, :math:`\ell=0`)
degenerate of the multi-mode scattering kernel :math:`R\circ\Lambda\circ
M`; a :class:`~orpheus.numerics.frame.FrameBase` manages the analogous
*stack* of dyads. The 0-ULP equivalence with the matvec arm is pinned by
``tests/sn/operators/test_fission_kernel_crosscheck.py`` (B.2), the
per-term closed-form correctness by
``tests/transport/test_reaction_rate_functional.py``.

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

apply + a working ``apply_transpose`` (``is_adjointable=True``). No
``solve``: the operator is **rank-1 in energy** per cell — its inverse does
not exist. The adjoint fission operator :math:`F^\dagger\,\phi^* =
\nu\Sigma_f \cdot (\chi \cdot \phi^*)` IS advertised (campaign #276): it is
the **dual dyad** :math:`|\nu\Sigma_f\rangle\langle\chi|` — the χ↔νΣf role
swap (emission becomes the adjoint weighting, production the adjoint
spectrum). It falls out of :attr:`kernel`'s tensor-product transpose by
swapping the rank-1 ``outer`` factors — no new kernel code; the transpose
lives on the :class:`~orpheus.numerics.operator.RankOneOperator` primitive.
See :meth:`apply_transpose`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import singledispatchmethod
from typing import TYPE_CHECKING, Any, overload

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    IdentityOperator,
    LinearOperator,
    TensorProductOperator,
    outer,
)

# Runtime imports for :func:`singledispatchmethod.register` — see
# ``scattering.py`` for the same pattern.  These types form a leaf in
# the SN dependency graph (they do not import fission.py).  The L2
# pure-Field :class:`AngularFlux` is re-aliased to ``AngularFlux`` to
# disambiguate from the legacy ``orpheus.sn.angular_flux.AngularFlux``
# which still rides on the operator-algebra path until D-H.1c.
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.reaction_rate_functional import ReactionRateFunctional
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
        ``{"apply", "apply_transpose"}`` — the rank-1 structure forbids a
        useful inverse (no ``solve``), but the dyad HAS a transpose: the
        adjoint fission :math:`F^\dagger = |\nu\Sigma_f\rangle\langle\chi|`
        (campaign #276), the χ↔νΣf role swap.
    """

    mat_xs: "MaterialXSField"


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

    @property
    def is_adjointable(self) -> bool:
        # F = outer(χ, ⟨νΣ_f|) is a rank-1 dyad; F^T is the free dyad-swap
        # (#112). is_invertible inherits base
        # False — a rank-1 production operator is singular.
        return True

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
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.full_field_space` instance
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

            \texttt{outer}\bigl(\chi,\;
                \mathrm{ReactionRateFunctional}(\nu\Sigma_f)\bigr)
            \;\&\; \mathrm{IdentityOperator}()

        which composes the rank-1 dyad :math:`F = |\chi\rangle\langle
        \nu\Sigma_f|` (Grand Report v3 §15 line 2008) as a type-visible
        separable form per §16A.10's ``B = G_patch \otimes K_omega
        \otimes K_g`` decomposition.  The first factor is the rank-1
        dyad — its reconstruction column is :math:`\chi` and its row
        co-vector is :attr:`production_rate` (the
        :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
        over :math:`\nu\Sigma_f`); the second factor advertises the
        spatial-axis broadcast.

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
            The 2-factor TP whose first element is the rank-1 dyad
            ``outer(self.chi, self.production_rate)`` (a
            :class:`RankOneOperator`) and whose second element is an
            :class:`IdentityOperator`.

        Notes
        -----
        The verification gate for Wave T step T.2 inspects this
        property — the type signature change from the legacy
        opaque-einsum body to a typed :class:`TensorProductOperator`
        is the deliverable.  See
        ``tests/sn/operators/test_fission_operator.py::TestRankOneTensorProductKernel``.
        """
        return outer(self.chi, self.production_rate) & IdentityOperator()

    @property
    def production_rate(self) -> ReactionRateFunctional:
        r"""The §5.6 production-rate co-vector — fission's rank-1 ROW-FACTOR.

        Returns the
        :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
        over this operator's production cross section :math:`\nu\Sigma_f` (read
        through the typed accessor
        :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.fission_production_field`).
        It contracts the group axis of a flux to the per-cell fission emission
        density

        .. math::

            p(\vec r) \;=\; \langle \nu\Sigma_f, \phi\rangle
                       \;=\; \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r) .

        This is the **row-factor of fission's rank-1 dyad**: :attr:`kernel` is
        literally ``outer(self.chi, self.production_rate) & IdentityOperator()``,
        so :meth:`apply` **routes through** this functional's
        :meth:`~orpheus.numerics.functional.InnerProductFunctional.evaluate` —
        the matvec's contraction IS this co-vector, not a parallel description
        of one (the procedural twin is dissolved). Together with :attr:`kernel`
        it makes :class:`FissionOperator` satisfy the §5.6
        :class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
        Protocol, and exposes the most physically central diagnostic in
        criticality — the per-cell fission source density — as a typed,
        inspectable :class:`~orpheus.numerics.functional.Functional`
        (``coding-elegance`` Pattern 3). The SAME type with the absorption cross
        section is the absorption rate :math:`\langle\Sigma_a,\phi\rangle`; the
        Rayleigh-quotient eigenvalue ``k = production / absorption`` is their
        ratio.

        Fission is the rank-1 (single-mode) degenerate of the multi-mode
        scattering kernel ``R∘Λ∘M``: the production-rate co-vector is the
        :math:`\ell=0` analysis row and :math:`\chi` the reconstruction column.
        Built fresh on each access (O(1) — one
        :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
        wrap) to honour the :attr:`mat_xs` read-through: a depletion /
        thermal-feedback update to :math:`\nu\Sigma_f` shows up immediately in
        the next :meth:`apply`.

        Correctness is pinned structurally against the closed-form
        :math:`k_\infty = \lambda_{\max}(A^{-1}F)` in
        ``tests/transport/test_reaction_rate_functional.py`` (production and
        absorption each, independently); the 0-ULP equivalence with the matvec
        arm is ``tests/sn/operators/test_fission_kernel_crosscheck.py`` (B.2).
        """
        return ReactionRateFunctional(self.mat_xs.fission_production_field)

    def apply_transpose(self, phi_star: np.ndarray) -> np.ndarray:
        r"""Adjoint fission :math:`F^\dagger\psi^* = \nu\Sigma_f\,(\chi\cdot\psi^*)`.

        The transpose of the rank-1 dyad
        :math:`F = |\chi\rangle\langle\nu\Sigma_f|` is the **dual dyad**
        :math:`F^\dagger = |\nu\Sigma_f\rangle\langle\chi|` — the χ↔νΣf **role
        swap**: the emission spectrum :math:`\chi` becomes the contracted row
        co-vector and the production cross section :math:`\nu\Sigma_f` becomes
        the reconstruction column,

        .. math::

            (F^\dagger\psi^*)_g(\vec r) = \nu\Sigma_{f,g}(\vec r)\,
                \sum_{g'} \chi_{g'}(\vec r)\,\psi^*_{g'}(\vec r).

        Physically: a high adjoint flux (importance) in the emission groups
        :math:`\chi` makes a cell a strong adjoint source weighted by its
        production :math:`\nu\Sigma_f` — the adjoint of "fission in :math:`g'`
        emits into the :math:`\chi` spectrum".

        This is the Euclidean transpose :math:`F^{T}`; the metric-correct
        Hilbert adjoint :math:`F^\dagger = G^{-1}F^{T}G` is the
        :attr:`~orpheus.numerics.operator.LinearOperator.H` wrapper's job.
        Routes through :attr:`kernel`'s
        :meth:`~orpheus.numerics.operator.TensorProductOperator.apply_transpose`
        — single source of truth: the SAME rank-1 ``outer`` primitive with its
        column and row swapped
        (:meth:`~orpheus.numerics.operator.RankOneOperator.apply_transpose`),
        the :class:`~orpheus.numerics.operator.IdentityOperator` spatial factor
        self-transposing.

        Bare-``np.ndarray`` surface — the K-eigenvalue / adjoint outer-iteration
        boundary, symmetric with the bare-ndarray :meth:`apply` arm. Accepts a
        :class:`~orpheus.numerics.field.Field`-like carrier (its ``.values`` is
        read); typed adjoint-flux carriers are threaded when the adjoint
        transport machinery lands (campaign #276).
        """
        arr = np.asarray(getattr(phi_star, "values", phi_star))
        return self.kernel.apply_transpose(arr)

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
          implicit-zero :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`
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
        (history-blind). The ``@overload`` static stubs name ``FullField``
        too (W-F), matching this runtime registration.

        Math: identical to the :class:`AngularFlux` branch above —
        reduce bulk angular → scalar via :math:`\phi = \sum_n w_n
        \psi_n`, compute iso fission source :math:`F\phi`, project to
        per-ordinate via :math:`/W`.  The output bulk is a pure-Field
        :class:`AngularFlux`; the output boundary is an **implicit
        zero** :class:`AngularBoundaryFlux` because the fission operator
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
        # Parse the arm's real contract per family (#289 discipline). The
        # ANGULAR arm consumes the PER-ORDINATE bulk (it needs
        # ``integrate_angular``) and is AngularFlux-only — deliberately
        # NARROWER than scattering's composite arm (which also accepts a
        # windowed HarmonicMomentFlux bulk). The SCALAR arm (#290 P4) is
        # the diffusion / CP composite: the bulk IS already the scalar
        # flux, so the reduction is the identity and the iso source needs
        # no per-ordinate projection. Both arms reuse the ScalarFlux
        # branch — single source of truth for the per-cell
        # production-rate × emission-spectrum contraction.
        bulk = psi.bulk
        if isinstance(bulk, AngularFlux):
            phi_scalar = bulk.integrate_angular()
            fission_iso: ScalarSourceSink = self.apply(phi_scalar)
            from orpheus.transport.source_sinks import (
                AngularSourceSink,
                AngularBoundarySourceSink,
            )
            # The composite's one mesh, read off the PARSED angular bulk —
            # ``AngularField.mesh`` carries the SNMesh declaration the widened
            # composite surfaces erase (#290 P2; same object by the composite's
            # mesh-identity invariant).
            mesh = bulk.mesh
            per_ord = AngularSourceSink.from_isotropic(
                fission_iso.values, mesh,
            )
            return FullField(
                # B.5.2: the operator output IS a source (Fψ rate density) — emit
                # the AngularSourceSink directly, not a re-wrap into AngularFlux.
                # #257 S8a: history-free (the matvec leaf is a base arrow; the
                # comonad lives on the driver, which reattaches the timed type when
                # this source is added to the timed rhs).
                bulk=per_ord,
                boundary=AngularBoundarySourceSink.zeros_on(mesh),
            )
        if isinstance(bulk, ScalarFlux):
            # Scalar composite arm (#290 P4): fission emission in iso
            # scalar magnitude on the (J⁺, J⁻)-trace composite — the
            # ScalarFlux branch supplies the bulk, the boundary is the
            # implicit-zero scalar source/sink (fission has no boundary
            # action; χ lives on cell-centred volumes).
            from orpheus.transport.source_sinks import ScalarBoundarySourceSink

            scalar_iso: ScalarSourceSink = self.apply(bulk)
            return FullField(
                bulk=scalar_iso,
                boundary=ScalarBoundarySourceSink.zeros_on(bulk.mesh),
            )
        raise TypeError(
            f"FissionOperator's composite arm requires an AngularFlux "
            f"bulk (the angular reduction φ = ∫ψ dΩ) or a ScalarFlux "
            f"bulk (the scalar composite, #290 P4); got "
            f"{type(bulk).__name__}."
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
        def apply(self, phi: FullField, /) -> "FullField": ...
        @overload
        def apply(self, phi: ScalarFlux, /) -> "ScalarSourceSink": ...
        @overload
        def apply(self, phi: np.ndarray, /) -> np.ndarray: ...
        def apply(self, phi: Any, /) -> Any: ...
    else:
        apply = _apply_impl

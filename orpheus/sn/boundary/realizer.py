r"""SN-specific BoundaryRealizer.

Dispatches by ``isinstance(law, ...)`` to the Wave-0 / Wave-1 primitives
that realise each :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
concrete as a single-arg :class:`~orpheus.numerics.operator.LinearOperator`
consumable by the SN sweep + Krylov path.

The isinstance dispatch keys on the Grand Report v3 canonical
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` concretes
(:class:`VacuumInflow`, :class:`ReflectiveBoundary`,
:class:`WhiteBoundary`, :class:`AlbedoBoundary`,
:class:`PeriodicBoundary`, :class:`PrescribedInflow`). The deprecated
``*BoundaryOperator`` aliases (Wave 5 ship-state) were retired in
Wave O step O.4a.1.

The realisation map (law concrete → Wave-0 / Wave-1 primitive)
================================================================

Since campaign phase **B3.2** a realized SN law is typed
:math:`\Gamma_+ \to \Gamma_-` — it consumes the **outflow** half-trace and
produces the **inflow** half-trace, exactly as the affine form
:math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q` says and as the diffusion arm
already did. The consumer composes ``ι₋ ∘ law ∘ γ₊``; nothing is computed and
then discarded. Realizing a law therefore needs the face's
:attr:`~orpheus.sn.mesh.method_space.SNMethodSpace.outflow_indices` (its
domain) alongside ``inflow_indices`` (its codomain) — face-orientation data
that a quadrature alone cannot supply, so ``SNMethodSpace.minimal`` no longer
suffices for the two laws below.

* :class:`~orpheus.geometry.boundary.vacuum.VacuumInflow` → the **zero map**
  :math:`\Gamma_+ \to \Gamma_-` (a :class:`~orpheus.numerics.operator.ZeroOperator`
  carrying both space hooks). Vacuum's whole content is :math:`R = 0`; with
  the domain narrowed there is nothing else to represent.

  Pre-B3.2 this was an :class:`IncomingOrdinateMaskTensor` — a full-face
  projector onto the OUTFLOW subspace whose preserved rows the consumer then
  discarded. Two campaign phases documented that survival as having "no
  consumer today"; the narrowing removes the question rather than answering
  it, because those rows are no longer in the operator's domain.
* :class:`~orpheus.geometry.boundary.reflective.ReflectiveBoundary(axis, albedo)` →
  ``albedo * PermutationOperator(local_perm)`` on the REDUCED ordinate axis,
  where ``local_perm = γ₊.to_local(reflection_index(axis)[inflow])`` — row
  :math:`j` of the image reads the mirror of the :math:`j`-th inflow ordinate,
  at that ordinate's position inside :math:`\Gamma_+`. (The remap must go
  through ``to_local``: on a slab the mirror REVERSES order, so a hand-written
  ``arange`` is wrong there — see the ``TraceRestrictionOperator`` docstring.)
  The ``albedo=1.0`` fast path returns the bare
  :class:`PermutationOperator` TP.
* :class:`~orpheus.geometry.boundary.white.WhiteBoundary(axis, outward_sign, albedo)` →
  ``albedo * (IsotropicEmissionOperator(...) @ PartialCurrentOperator(...))``
  (with the ``albedo=1.0`` fast path), narrowed at **B3.4a** to contract over
  :math:`\Gamma_+` and re-emit on :math:`\Gamma_-`. Two things dissolved with
  the codomain: the operator's private ``> 0.0`` outflow test (it now reads the
  one face-name → signed-projection primitive against ``TANGENTIAL_EPS``, so it
  can no longer disagree with the trace space), and the full-face broadcast that
  wrote outflow rows nobody read. The law's declared ``axis``/``outward_sign``
  is cross-checked against the installation face's :math:`\Gamma_+` — see
  :func:`_checked_angular_average`.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundary(α, closure)` →
  the SAME body its closure names, at **B3.4b**: a
  :class:`~orpheus.geometry.boundary.SpecularReturn` routes to
  :func:`_specular_kernel` (shared with reflective) and an
  :class:`~orpheus.geometry.boundary.IsotropicReturn` to
  :func:`_checked_angular_average` (shared with white). So
  ``AlbedoBoundary(α, SpecularReturn(a)) ≡ ReflectiveBoundary(a, α)`` and
  ``AlbedoBoundary(α, IsotropicReturn(a, s)) ≡ WhiteBoundary(a, s, α)`` as
  matrices, by executing one construction rather than by two transcriptions
  agreeing. The laws still assert different physics — a wall's constitutive
  return versus a symmetry of the domain — which is the user's 2026-08-01
  ruling: the specular pairing belongs to :math:`R`, and
  ``AlbedoBoundary.geometry_map`` is ``SelfPairedDeck.identity()``
  unconditionally.
* :class:`~orpheus.geometry.boundary.albedo.AlbedoBoundary(α)` with **no**
  closure → **REFUSED** (:class:`~orpheus.geometry.boundary.BoundaryError`).
  Its response :math:`R = \alpha\,I` is an endomorphism of :math:`\Gamma_+` and
  its :math:`G` supplies no crossing, so on an ANGULAR trace nothing says which
  outgoing direction feeds which incoming one; composing it anyway pairs them
  by array position. Pre-B3.4b that is exactly what happened — the arm returned
  ``ZeroOperator`` / ``IdentityOperator`` / ``α·(I & I)``, full-face
  endomorphisms which the ``ι₋ ∘ law ∘ γ₊`` composite then read positionally.
  The refusal is of the **incomplete spelling**, not of the law: a SCALAR
  method needs no closure (:math:`J^- = \alpha J^+` has one degree of freedom),
  so the diffusion realizer takes the same object unchanged and
  ``BC("albedo", albedo=…)`` is untouched. This is the method-realizability
  taxonomy's **angular-resolution** axis biting for the first time.
* :class:`~orpheus.geometry.boundary.periodic.PeriodicBoundary` →
  :class:`~orpheus.numerics.operator.IdentityOperator` (the crossing it
  carries is the PARTNER-face channel, not an action on the trace).
* :class:`~orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow(source)` →
  :class:`~orpheus.numerics.operator.ZeroOperator` — the same zero morphism
  vacuum returns, because the law is affine and this tier realizes its LINEAR
  factor, which is :math:`0`. The source travels the boundary-source channel
  (:ref:`bc-affine-source-channel`). The rank-0 affine BC (Wave 7 / C7.5);
  until **P3** it realized to an ``IncomingSourceOperator`` whose apply
  ignored its input and returned :math:`q`. **B3.4a** retired the inflow MASK
  (#52 / ERR-047) along with the codomain it corrected: the rows it zeroed are
  no longer emitted on, so :math:`q \in \Gamma_-` holds by typing rather than by
  erasure. Its unmasked companion branch retired with it — post-B3.2 a method
  space without face data cannot supply :math:`\gamma_+` either, so that
  fallback was already unreachable here.
* :class:`~orpheus.geometry.boundary.zero_flux.ZeroFluxBoundary` →
  **REFUSED** (:class:`~orpheus.geometry.boundary.BoundaryError`):
  the Dirichlet idealization :math:`\phi_\Gamma = 0` is the
  albedo-family member :math:`\mathcal{A} = -1`, and a negative
  angular inflow is unrepresentable in a non-negative :math:`\psi`.
  It realizes only under the diffusion realizer (#290 P3).

Realization is the certification point (#52): every dispatch first
fires the law's :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_realizable`
aggregate (the §16A.12 universal invariants + each law's specific
checks — reflective's measure/involution/inflow→outflow table trio,
white/albedo's sub-Markov bound, the universal ERR-047 source-trace
certification), and the vacuum arm cross-checks the claimed
``inflow_indices`` against the face-name geometry (ERR-041). The
error catalog's recurring lesson operationalized: contracts are
checked where the law meets the trace, not discovered by downstream
balance defects.

Wave 11 removed the ``MixedBoundaryOperator`` composer from the
dispatch table. Rank-N compositions are now expressed via Wave-0
``OperatorSum``-algebra over already-realised leaves; the walker for
``BoundaryTraceLaw``-rooted expression trees is the method-blind
:func:`~orpheus.geometry.boundary.realize_recursively` (moved to
``geometry/boundary/`` at #290 P7b — pass
``realizer=SNBoundaryRealizer()`` at the leaves).

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 5 — C5.3 (functional
  realizer brief), §16A.5 (vacuum-trace semantic correction), Wave 11
  (mixed-BC removal + tree-walking realisation).
* Grand Report v3 §16A.3 lines 2841–2860 (realizer-as-third-layer
  motivation), §16A.10 lines 3160-3175 (vacuum-trace correct
  representation).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    BoundaryTraceLaw,
    LambertianReemission,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    ScalarResponse,
    SpecularReemission,
    VacuumAppliedToOutgoingTraceError,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
    stamp_boundary_role,
)
from orpheus.numerics.face_layout import AXIS_NAMES, face_name
from orpheus.numerics.spaces.angular_trace_space import (
    TANGENTIAL_EPS,
    build_omega_dot_n,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
    PermutationOperator,
    TraceRestrictionOperator,
    ZeroOperator,
    checked_space_extent,
)
from .angular import (
    IsotropicEmissionOperator,
    PartialCurrentOperator,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from orpheus.geometry.boundary import PairedDeck
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.mesh.method_space import SNMethodSpace


__all__ = ["SNBoundaryRealizer"]


def _zero_rows(n_rows: int) -> "Callable[[object], np.ndarray]":
    r"""A :class:`ZeroOperator` space hook emitting ``n_rows`` ordinates.

    The zero map between two DIFFERENT spaces must emit the zero of the space
    it lands in, not an echo of the one it came from — see
    :class:`~orpheus.numerics.operator.ZeroOperator`'s own note. Trailing axes
    (group, spatial) are carried through unchanged.

    The parameter is typed ``object`` because a hook CONSUMED by the operator
    is contravariant in its input: it must accept whatever the operator hands
    it. This mirrors ``_ray_source_zero`` in :mod:`orpheus.sn.solver`, the
    codebase's existing space-hook spelling.
    """
    def build(x: object) -> "np.ndarray":
        arr = np.asarray(x)
        return np.zeros((n_rows,) + arr.shape[1:], dtype=arr.dtype)

    return build


def _outflow_restriction(
    method_space: "SNMethodSpace", law_key: str,
) -> "TraceRestrictionOperator":
    r"""The face's :math:`\gamma_+` — a narrowed law's DOMAIN.

    Since B3.2 an SN boundary law is typed :math:`\Gamma_+ \to \Gamma_-`, so
    realizing one requires knowing which ordinates are outflow AT THIS FACE.
    That is face-orientation data, not quadrature data: a law realized without
    a face has no :math:`\Gamma_+` to restrict to, and the demand is
    structural rather than an implementation detail. (Pre-B3.2 the specular
    branch needed only ``quadrature.reflection_index``, which is why a
    quadrature-only ``SNMethodSpace.minimal`` sufficed for it.)

    Built here rather than fetched from the trace's cache because realization
    runs **once per face at mesh construction**, not per sweep — the cached
    accessor exists for the per-sweep consumer. Reading the method space's own
    field also keeps a hand-built space workable, which is the reason
    ``inflow_indices`` is a field too.

    Bound :math:`\Gamma(f) \to \Gamma_+(f)` when the method space carries a
    trace (G6.3 step 5, #330). Building here rather than fetching the cached
    :math:`\gamma_+` is a construction-cost decision, not a typing one — so
    the local build owes the same binding the cached accessor already has, or
    the realizer would hand its own consumer an untyped domain while typing
    that consumer's output.
    """
    if method_space.outflow_indices is None:
        raise BoundaryError(
            f"SNBoundaryRealizer cannot realize {law_key!r} without "
            f"outflow_indices: since B3.2 a boundary law is typed Γ₊ → Γ₋, "
            f"and its DOMAIN is the outflow half-trace of a particular face — "
            f"which a quadrature alone cannot name. Construct via "
            f"SNMethodSpace.for_face(quadrature=..., face=..., trace=...), "
            f"which populates both half-traces, or pass outflow_indices "
            f"explicitly alongside inflow_indices.",
            law=law_key,
        )
    trace, face = method_space.trace, method_space.face
    gamma, gamma_plus = (
        (trace.face_space(face), trace.outflow_space(face))
        if trace is not None and face is not None
        else (None, None)
    )
    return TraceRestrictionOperator(
        np.sort(np.asarray(method_space.outflow_indices, dtype=np.intp)),
        n_total=method_space.quadrature.N,
        axis=0,
        domain=gamma,
        codomain=gamma_plus,
    )


def _narrowed_zero_operator(
    method_space: "SNMethodSpace",
    gamma_out: "TraceRestrictionOperator",
    *,
    law_key: str,
) -> "ZeroOperator":
    r"""The zero map :math:`\Gamma_+ \to \Gamma_-`, with both spaces named.

    A zero map between two DIFFERENT spaces must emit the zero of the space it
    lands in, not an echo of the one it came from — see
    :class:`~orpheus.numerics.operator.ZeroOperator`'s own note. The forward
    emits the zero of :math:`\Gamma_-`, the transpose the zero of
    :math:`\Gamma_+`. Relying on the endomorphic ``0.0 * x`` echo would be
    wrong in principle and merely lucky in practice: ``|Γ₊| == |Γ₋|`` on every
    reachable fixture, so the shapes coincide — an accident, not a contract.

    THREE callers, one body: vacuum (whose response IS structurally zero), any
    re-emission law at :math:`\alpha = 0` (whose response *evaluates* to zero —
    a perfectly absorbing wall is a vacuum, and says so with the same object),
    and — since **P3** — prescribed inflow, whose LINEAR factor is zero and
    whose entire content is the affine :math:`q`
    (:ref:`bc-affine-source-channel`). That third caller is the one that makes
    the shared body a statement rather than a coincidence: vacuum and prescribed
    are the SAME operator, and differ only in a term this tier does not carry.

    Both half-traces are REQUIRED, and the demand is the whole point: a zero map
    that does not know its codomain degenerates to exactly the endomorphic echo
    described above.

    Since **G6.3 step 6** the two spaces are NAMED as well as sized. :math:`\Gamma_+`
    is read off ``gamma_out``'s own codomain rather than fetched from the trace a
    second time — ``_outflow_restriction`` bound it at step 5, and re-deriving it
    here would be a twin source for one fact. Each space is checked against the
    row count the matching hook was sized from: the space and the length describe
    the same fact, and a disagreement is invisible at apply-time because the zero
    arrays broadcast either way.
    """
    if method_space.inflow_indices is None:
        raise BoundaryError(
            f"SNBoundaryRealizer cannot build the zero map for {law_key!r} "
            f"without inflow_indices: the map is Γ₊ → Γ₋ and must emit the "
            f"zero of Γ₋, which it cannot size. Construct via "
            f"SNMethodSpace.for_face(quadrature=..., face=..., trace=...).",
            law=law_key,
        )
    n_inflow = int(method_space.inflow_indices.size)
    trace, face = method_space.trace, method_space.face
    gamma_minus = (
        trace.inflow_space(face)
        if trace is not None and face is not None
        else None
    )
    return ZeroOperator(
        # ``reportArgumentType``: ``ZeroOperator`` is generic over ``Vector``,
        # and MEASURED (2026-07-31, pyright + numpy stubs) ``np.ndarray`` does
        # not satisfy that protocol statically — ``Self``-typed ``__add__``
        # will not bind against numpy's overloads. The gap is upstream, not at
        # this call site: every ndarray-level operator in ``numerics.operator``
        # sidesteps it by being declared UNPARAMETERIZED
        # (``PermutationOperator(LinearOperator)``), an option a generic's own
        # hooks do not have. Runtime conformance is real; only the static bind
        # fails.
        codomain_zero=_zero_rows(n_inflow),  # type: ignore[reportArgumentType]
        transpose_zero=_zero_rows(  # type: ignore[reportArgumentType]
            gamma_out.n_restricted
        ),
        domain=checked_space_extent(
            gamma_out.codomain, gamma_out.n_restricted, axis=0,
            owner="_narrowed_zero_operator", role="domain",
        ),
        codomain=checked_space_extent(
            gamma_minus, n_inflow, axis=0,
            owner="_narrowed_zero_operator", role="codomain",
        ),
    )


def _attenuated_kernel_operator(
    k_omega: "LinearOperator",
    alpha: float,
    *,
    method_space: "SNMethodSpace",
    gamma_out: "TraceRestrictionOperator",
    law_key: str,
) -> "LinearOperator":
    r"""``α · (K_ω ⊗ K_g)``, stamped :attr:`BlockRole.BOUNDARY`.

    The §16A.10 BC decomposition is ``B = G_patch ⊗ K_omega ⊗ K_g``; every
    narrowed leaf this realizer builds is that shape with the group factor an
    identity and the amplitude out front.

    One body, because the assembly is genuinely one thing: before **B3.4b**
    this pattern was written once in the reflective arm and once in the white
    arm, and B3.4b's two albedo closures would have made four transcriptions of
    it. All four kernel routes end here instead.

    The three amplitude regimes, and why the ends are not merely optimizations
    ------------------------------------------------------------------------

    * :math:`\alpha = 1` returns the bare tensor product, so the fold reduces
      to the angular kernel's own action with no ``ScaledOperator`` wrapper.
      A convenience.
    * :math:`\alpha = 0` returns :func:`_narrowed_zero_operator`, because a
      surface that returns nothing IS a vacuum — the same object, honestly
      typed, with a working transpose. **Not** a convenience: the general path
      cannot express it, since
      :class:`~orpheus.numerics.operator.ScaledOperator` refuses a zero scalar
      as degenerate. Before B3.4b that refusal was reachable —
      ``ReflectiveBoundary(axis, 0.0)`` and ``WhiteBoundary(..., 0.0)`` are
      legal laws (:math:`\alpha = 0` satisfies every invariant, including
      sub-Markov) and both died in the numerics layer with a message about
      operator degeneracy rather than realizing the boundary they describe.
      Folding the four routes into one body is what made that visible, and one
      answer fixes all four.
    * :math:`0 < \alpha < 1` takes the scaled path.

    .. note::

       At :math:`\alpha = 0` the caller has already built ``k_omega`` and this
       function discards it. That is deliberate, and it is not the
       compute-then-discard pattern this campaign removes: what is kept is the
       **construction's validation**, not its value. Building the specular
       kernel runs ``γ₊.to_local``, which refuses a mirror that sends an inflow
       ordinate outside :math:`\Gamma_+` (ERR-045 at realization). Skipping it
       would make *realizability itself* depend on the amplitude — a law that
       refuses at :math:`\alpha = 0.001` and is accepted at :math:`\alpha = 0`
       — and a declared-but-unrealizable pairing is a nonsense declaration at
       every amplitude.
    """
    if alpha == 0.0:
        return stamp_boundary_role(
            _narrowed_zero_operator(method_space, gamma_out, law_key=law_key)
        )
    base = k_omega & IdentityOperator()  # 2-factor TensorProductOperator
    if alpha == 1.0:
        return stamp_boundary_role(base)
    return stamp_boundary_role(float(alpha) * base)


def _specular_kernel(
    quadrature: "Quadrature",
    method_space: "SNMethodSpace",
    gamma_out: "TraceRestrictionOperator",
    *,
    axis: str,
    law_key: str,
) -> "PermutationOperator":
    r"""The specular pairing on :math:`\Gamma_+ \to \Gamma_-`, as a permutation
    of the REDUCED ordinate axis.

    The mirror sends the inflow ordinate ``inflow[j]`` to the outflow ordinate
    ``perm[inflow[j]]`` (the ERR-045 invariant the law layer certifies at
    realization: a reflection maps inflow to OUTFLOW, never to itself). So the
    narrowed operator's row ``j`` reads the position of that outflow ordinate
    INSIDE :math:`\Gamma_+` — which is exactly what ``γ₊.to_local`` computes,
    and what a hand-written ``arange`` would get wrong: on a slab the mirror
    REVERSES order, so ``perm[inflow] = [3, 2]`` maps to local ``[1, 0]``.

    **Two laws reach this body** (B3.4b), carrying the pairing in different
    tiers: :class:`~orpheus.geometry.boundary.ReflectiveBoundary` in :math:`G`
    (a symmetry of the domain) and
    :class:`~orpheus.geometry.boundary.AlbedoBoundary` with a
    :class:`~orpheus.geometry.boundary.SpecularReturn` closure in :math:`R` (a
    polished wall). They assert different physics and realize to the same
    matrix — so they share this construction rather than agreeing by
    transcription. ``law_key`` names the caller in the errors only.

    The deck transformation as a length-1 chain
    -------------------------------------------

    Returned **bound** to :math:`\Gamma_+(f) \to \Gamma_-(f)` (G6.3 step 5,
    #330), which makes it the degenerate case of the same structure the
    diffuse arm builds in :func:`_checked_angular_average`: a boundary law is
    a chain from outflow to inflow, and a measure-preserving bijection has
    nothing to factor, so its chain has ONE link. There is no separate
    "atomic" code path — only a shorter chain.

    ⭐ **Binding is also what retired the involution flag.** `[M]` the
    narrowed local permutation satisfies ``perm[perm] == arange`` on
    ``gauss_legendre(4/8)``, ``product(4,4)`` and ``level_symmetric(6)``, and
    NOT on ``lebedev(17)`` — one physical law, an answer tracking the
    quadrature's local index ordering rather than the mirror. Bound, the
    question stops being answerable at all: :math:`P \circ P` is not an
    expression when the ends are different spaces, so
    :class:`~orpheus.numerics.operator.PermutationOperator` no longer stores
    a ``bool`` about it — ``P @ P`` raises, which is the same claim delivered
    by the algebra. The involution that IS real lives one tier up, on the
    full-space table
    (:meth:`~orpheus.numerics.quadrature.Quadrature.reflection_index`), where
    domain and codomain coincide and ERR-044 guards it.
    """
    inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
    perm = quadrature.reflection_index(axis)
    if inflow.size != gamma_out.n_restricted:
        raise BoundaryError(
            f"SNBoundaryRealizer cannot narrow a specular pairing at "
            f"face {method_space.face!r}: |Γ₋| = {inflow.size} but "
            f"|Γ₊| = {gamma_out.n_restricted}. A specular mirror is a "
            f"BIJECTION between the two half-traces, so a face where "
            f"they differ in size has no specular realization — which "
            f"means the quadrature's tangential band has swallowed "
            f"ordinates asymmetrically at this face.",
            law=law_key,
        )
    # ``to_local`` raises if the mirror sent an inflow ordinate anywhere but
    # Γ₊ — the ERR-045 violation, caught here as a crossed-index-set error
    # rather than as silent wrong physics. It stays the SINGLE authority on
    # that index question (no second membership test here to drift from it);
    # what this arm adds is the semantic diagnosis, so the caller gets an
    # attributed BoundaryError naming the two possible causes rather than a
    # raw index complaint from a helper it never called. Symmetric with the
    # diffuse arm's orientation cross-check — before B3.4b these two
    # disagreed, because a law had no way to declare an axis that could
    # disagree with its installation face until the closure gave it one.
    try:
        local_perm = gamma_out.to_local(perm[inflow])
    except ValueError as exc:
        raise BoundaryError(
            f"A specular pairing about axis {axis!r} does not map Γ₋ into Γ₊ "
            f"at face {method_space.face!r}. Two causes: the declared axis "
            f"disagrees with the face the law is installed on — a mirror "
            f"about 'y' on an x-face relabels WITHIN each half-trace instead "
            f"of exchanging them, which is not a boundary law at all — or "
            f"the quadrature's reflection table violates ERR-045 (an inflow "
            f"ordinate whose partner is not outflow). The law-level "
            f"certification catches the second at assert_realizable, so at "
            f"this point the first is the likely one. Underlying: {exc}",
            law=law_key,
        ) from exc

    # Bind the single link to the Γ ladder when the method space carries a
    # trace — the canonical `SNMethodSpace.for_face` path always does. A
    # hand-built space may not, and binding stays OPTIONAL until the tree-wide
    # mandate (#330): an unbound permutation gathers the same rows, it just
    # forfeits the composability check and the metric-aware `.H`.
    trace, face = method_space.trace, method_space.face
    gamma_plus = gamma_minus = None
    if trace is not None and face is not None:
        gamma_plus = trace.outflow_space(face)
        gamma_minus = trace.inflow_space(face)
    return PermutationOperator(
        local_perm, axis=0, domain=gamma_plus, codomain=gamma_minus,
    )


def _assert_wrap_identification(
    wrap: "PairedDeck",
    method_space: "SNMethodSpace",
) -> str:
    r"""Certify the quotient reading :math:`\Gamma_+(f') \equiv \Gamma_-(f)`,
    and return the partner face.

    The user's **B3.4c** ruling: build the partner-face channel now, and let
    the quotient reading be *asserted at realization* rather than baked into
    the mesh topology — so this identification is a guard, not a restructure.
    (The alternative was to make the mesh carry a face-partner map, which
    commits every mesh to a topology decision that belongs to a law.)

    It is a geometric theorem, not a coincidence: the two faces' outward
    normals are opposite, so a direction OUTGOING at :math:`f'` is INCOMING at
    :math:`f`. `[M]` measured 2026-08-01 as an exact set equality on
    ``gauss_legendre(8)``, ``product(2,4)``, ``level_symmetric(6)`` and
    ``lebedev(17)``, on both axis pairs. It is asserted anyway, for the reason
    every guard in this campaign is asserted: the theorem is about the
    *continuous* normals, and what the realization needs is a statement about
    THIS quadrature's classified index sets, which a tangential band or a
    hand-built method space can break.

    The comparison is against index SETS, never sizes — ``|Γ₊| == |Γ₋|`` on
    every quadrature in the tree, so a size check is Mode-12 blind (the same
    trap the diffuse arm's cross-check documents).

    Both encodings are sourced independently: the installation face's
    :math:`\Gamma_-` comes from the method space (built by the trace space),
    the partner's :math:`\Gamma_+` from the one face-name → signed-projection
    primitive applied to the name :meth:`~orpheus.geometry.boundary.PairedDeck.domain_face` derives. That
    is the ERR-041 discipline — two encodings of one orientation, compared.

    The face demand is checked FIRST, ahead of the shared
    :func:`_outflow_restriction`, so a faceless method space hears the reason
    that is true of periodic and of no other law — *your domain is a different
    face* — rather than the generic missing-``outflow_indices`` complaint every
    narrowed law shares. Same data missing, but only one message tells the
    caller what periodic actually needs.
    """
    face = method_space.face
    if face is None:
        raise BoundaryError(
            "SNBoundaryRealizer cannot realize PeriodicBoundary without a "
            "face: this is the one law whose DOMAIN is a DIFFERENT face's Γ₊ "
            "(γ₋ψ|_f = γ₊ψ|_f'), so the partner cannot be named without "
            "knowing where the law is installed. Construct via "
            "SNMethodSpace.for_face(quadrature=..., face=..., trace=...).",
            law="periodic",
        )
    if method_space.inflow_indices is None:
        raise BoundaryError(
            f"SNBoundaryRealizer cannot certify the periodic identification "
            f"at face {face!r} without inflow_indices: the claim being "
            f"checked is Γ₊(partner) == Γ₋(this face), and Γ₋ is the half "
            f"the method space carries.",
            law="periodic",
        )
    partner = wrap.domain_face(face)          # raises on a mis-declared axis
    gamma_out = _outflow_restriction(method_space, "periodic")
    omega_dot_n = build_omega_dot_n(method_space.quadrature, (partner,))[0]
    partner_outflow = np.flatnonzero(omega_dot_n > +TANGENTIAL_EPS)
    face_inflow = np.sort(
        np.asarray(method_space.inflow_indices, dtype=np.intp)
    )
    if not np.array_equal(partner_outflow, face_inflow):
        raise BoundaryError(
            f"A periodic law on face {face!r} identifies it with partner "
            f"{partner!r}, but Γ₊({partner}) = {partner_outflow.tolist()} is "
            f"not Γ₋({face}) = {face_inflow.tolist()}. The wrap is realized as "
            f"the identity on the ordinate index, which is sound only when the "
            f"partner's outgoing directions ARE this face's incoming ones — "
            f"the opposite-normals theorem. A mismatch means the two faces do "
            f"not carry opposite normals as classified on this quadrature, so "
            f"the pair is not a translation quotient and needs an explicit "
            f"gluing map (issue #178).",
            law="periodic",
        )
    if gamma_out.n_restricted != face_inflow.size:
        raise BoundaryError(
            f"A periodic law on face {face!r} has |Γ₊({face})| = "
            f"{gamma_out.n_restricted} but |Γ₋({face})| = {face_inflow.size}. "
            f"A translation quotient is a BIJECTION of half-traces, so a face "
            f"where they differ in size has no periodic realization — the "
            f"quadrature's tangential band has swallowed ordinates "
            f"asymmetrically at this face.",
            law="periodic",
        )
    return partner


def _checked_angular_average(
    quadrature: "Quadrature",
    method_space: "SNMethodSpace",
    gamma_out: "TraceRestrictionOperator",
    *,
    axis: str,
    outward_sign: int,
    law_key: str,
) -> "LinearOperator":
    r"""The Lambertian kernel, with the DECLARED orientation cross-checked
    against the face it is being installed on (the ERR-041 pattern).

    A diffuse law declares its own ``axis`` / ``outward_sign`` — on the law
    itself for :class:`~orpheus.geometry.boundary.WhiteBoundary`, on the
    closure for
    :class:`~orpheus.geometry.boundary.AlbedoBoundary` +
    :class:`~orpheus.geometry.boundary.IsotropicReturn` — while the method
    space independently names the face. Those are two encodings of ONE
    orientation, and until B3.4a nothing compared them: a white law declared
    for ``x`` and installed on ``ymin`` averaged over the wrong hemisphere and
    reported nothing. The check is the same shape as the vacuum arm's
    (``ERR-041``: swapped inflow/outflow face annotation), against index SETS
    rather than sizes, because ``|Γ₊| == |Γ₋|`` on every quadrature in the tree
    makes a size comparison Mode-12 blind.

    On the canonical ``SNMesh.realize_boundary_law`` path both encodings derive
    from the same face label, so the guard is green by construction; it bites
    on hand-built method spaces and on a mis-declared law.

    Returns the response as the **two-link chain**
    :math:`\Gamma_+(f) \xrightarrow{C} S(f) \xrightarrow{B} \Gamma_-(f)`
    (G6.3 step 3, #330) rather than the welded ``AngularAverageOperator`` it
    replaced: the intermediate is then the named outgoing partial current
    :math:`J^+` instead of an anonymous local, and the composite carries a
    transpose the welded form withheld. `[M]` bit-identical to the predecessor
    on every shipped quadrature.
    """
    declared_face = face_name(AXIS_NAMES.index(axis), outward_sign)
    omega_dot_n = build_omega_dot_n(quadrature, (declared_face,))[0]
    declared_outflow = np.flatnonzero(omega_dot_n > +TANGENTIAL_EPS)
    if not np.array_equal(declared_outflow, gamma_out.indices):
        raise BoundaryError(
            f"A diffuse re-emission law declares axis={axis!r}, "
            f"outward_sign={outward_sign:+d} — i.e. the face "
            f"{declared_face!r}, whose Γ₊ is {declared_outflow.tolist()} — "
            f"but it is being installed where Γ₊ is "
            f"{np.asarray(gamma_out.indices).tolist()} "
            f"(face={method_space.face!r}). The Lambertian averages over "
            f"the hemisphere its OWN orientation names, so a mismatch "
            f"averages the wrong ordinates silently.",
            law=law_key,
        )
    inflow = np.flatnonzero(omega_dot_n < -TANGENTIAL_EPS)
    if inflow.size == 0:
        raise BoundaryError(
            f"A diffuse re-emission law on face {declared_face!r} has no "
            f"inflow ordinates — the re-emitted flux has no Γ₋ to land on. "
            f"The quadrature is degenerate for this face.",
            law=law_key,
        )
    cos_w = np.asarray(quadrature.weights, dtype=float)[declared_outflow] * (
        omega_dot_n[declared_outflow]
    )

    # Bind the chain to the Γ ladder when the method space carries a trace —
    # the canonical `SNMethodSpace.for_face` path always does. A hand-built
    # space may not, and binding stays OPTIONAL until the tree-wide mandate
    # (#330): an unbound chain still computes the same numbers, it just
    # forfeits the composability check and the metric-aware `.H`.
    trace, face = method_space.trace, method_space.face
    gamma_plus = current = gamma_minus = None
    if trace is not None and face is not None:
        gamma_plus = trace.outflow_space(face)
        current = trace.current_space(face)
        gamma_minus = trace.inflow_space(face)

    contraction = PartialCurrentOperator(
        cos_w, domain=gamma_plus, codomain=current,
    )
    emission = IsotropicEmissionOperator(
        float(cos_w.sum()), int(inflow.size),
        domain=current, codomain=gamma_minus,
    )
    return emission @ contraction


class SNBoundaryRealizer:
    r"""Functional realizer for SN BCs.

    Dispatches by ``isinstance(law, ...)`` to map each
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` concrete
    onto a Wave-0 / Wave-1 primitive composed with scalar
    amplitudes (per :ref:`bc-tensor-decomposition` and the
    realisation map in this module's docstring).

    The realizer instance carries no state; ``method_name`` is the
    only attribute. The :meth:`realize` method is stateless — it
    can be called freely from any context.
    """

    method_name: str = "SN"

    def realize(
        self,
        law: "BoundaryTraceLaw",
        method_space: SNMethodSpace,
    ) -> LinearOperator:
        """Realize ``law`` for SN as a 1-arg :class:`LinearOperator`.

        The concrete return is always a :class:`LinearOperator`
        subclass (a generic numerics primitive —
        :class:`TensorProductOperator` / :class:`ScaledOperator` /
        :class:`~orpheus.numerics.operator.ZeroOperator` / …); the narrower mixin type
        (vs the bare :class:`LinearOperator` Protocol) lets the
        :func:`~orpheus.geometry.boundary.stamp_boundary_role`
        ``block_role`` instance-stamp type-check and keeps the
        realized-law primitives assignable nominally rather than
        through the (invariant) Protocol.
        """
        if isinstance(law, ZeroFluxBoundary):
            # ── REFUSAL AXIS: state-cone / sign ──────────────────────────
            # NOT an angular-resolution refusal and NOT a topological one —
            # see `bc-method-realizability` in orpheus.geometry.boundary. SN resolves
            # than diffusion does and couples faces just fine; what it cannot
            # do is hold this VALUE. The diffusion trace is a signed partial
            # current, so 𝒜 = −1 lives in its state space; the SN trace is an
            # angular flux with ψ ≥ 0, and a negative inflow is outside that
            # cone. The same law is realizable there and unrealizable here for
            # a reason that has nothing to do with either method's fidelity.
            #
            # #290 P3 (user ruling 3): φ_Γ = 0 is a diffusion-level
            # Dirichlet idealization — albedo 𝒜 = −1 in the
            # partial-current basis. A transport BC prescribes the
            # inflow trace from the outflow trace, and a NEGATIVE
            # angular inflow is unrepresentable in a non-negative ψ,
            # so SN refuses rather than realizing something else.
            raise BoundaryError(
                "SNBoundaryRealizer cannot realize ZeroFluxBoundary: "
                "the zero-flux Dirichlet condition (φ_Γ = 0) is the "
                "albedo-family member 𝒜 = −1 in the partial-current "
                "basis, and a negative angular inflow has no transport "
                "realization (ψ ≥ 0 ⟹ J± ≥ 0). Use VacuumInflow for "
                "the physical zero-incoming law (J⁻ = 0), or realize "
                "zero-flux with the diffusion realizer.",
                law="zero_flux",
            )

        quad = method_space.quadrature

        # #52 — the §16A.12 invariants fire HERE, at the seam where the
        # law meets the actual quadrature + trace data (the error
        # catalog's recurring lesson: check at realization, not only by
        # downstream balance). One aggregate call; each law extends the
        # base template with its own invariants (reflective's three
        # table checks, white/albedo's sub-Markov bound, the universal
        # ERR-047 source-trace certification). Certification is a LAW
        # contract: a non-law object (and a LawSum/LawScaled tree,
        # which composes laws but is not one) falls through to the
        # loud dispatch-failure raise below instead.
        if isinstance(law, BoundaryTraceLaw):
            law.assert_realizable(
                quad, inflow_indices=method_space.inflow_indices
            )

        if isinstance(law, VacuumInflow):
            if method_space.inflow_indices is None:
                raise BoundaryError(
                    "SNBoundaryRealizer cannot realize "
                    "VacuumInflow without inflow_indices in "
                    "method_space. Wave 8's SNMesh wiring populates "
                    "this from "
                    "AngularTraceSpace.inflow_indices_for_face. For "
                    "now, supply inflow_indices explicitly via "
                    "SNMethodSpace(quadrature=quad, face=..., "
                    "inflow_indices=...).",
                    law="vacuum",
                )
            # #52 / ERR-041 — trace-orientation guard. The claimed
            # inflow_indices and the face NAME are two independent
            # encodings of the same orientation (the annotation-swap
            # bug class desynchronizes them: indices built from the
            # OUTGOING set of this face). Cross-check against the
            # signed projection the face name alone implies via the
            # single face→normal primitive. On the canonical
            # SNMesh.realize_boundary_law path both encodings derive
            # from the trace space, so the guard is green by
            # construction; it bites on hand-built method spaces —
            # exactly where the annotation swap lives. A method space
            # with no face carries no independent orientation truth,
            # so the guard cannot fire there (documented escape:
            # quadrature-only spaces with explicit indices are the
            # caller's responsibility).
            if method_space.face is not None:
                omega_dot_n = build_omega_dot_n(
                    quad, (method_space.face,)
                )[0]
                claimed = np.asarray(
                    method_space.inflow_indices, dtype=np.intp
                )
                outgoing = omega_dot_n[claimed] > +TANGENTIAL_EPS
                if np.any(outgoing):
                    offending = claimed[outgoing]
                    raise VacuumAppliedToOutgoingTraceError(
                        f"VacuumInflow at face "
                        f"{method_space.face!r}: inflow_indices "
                        f"{offending.tolist()} are OUTGOING ordinates "
                        f"(Ω·n̂ > 0) at this face — the vacuum law "
                        f"prescribes γ₋ψ = 0 on the INCOMING trace "
                        f"only; zeroing outflow slots discards the "
                        f"flux the sweep just computed (ERR-041: "
                        f"swapped inflow/outflow face annotation).",
                        law="vacuum",
                    )
            # B3.2 — the domain is Γ₊, so vacuum is the honest ZERO MAP.
            #
            # Vacuum's whole content is ``R = 0``; its ``G`` is the identity
            # deck element (it fixes no geometry). With the law typed
            # ``Γ₊ → Γ₋`` there is nothing left to represent but the zero map
            # between those two spaces — and the question the old realization
            # raised, "why does a projector preserve rows nobody reads?",
            # dissolves rather than needing an answer: those rows are no
            # longer in the operator's domain at all.
            #
            # Pre-B3.2 this returned ``IncomingOrdinateMaskTensor(...) &
            # IdentityOperator()`` — a full-face projector onto the OUTFLOW
            # subspace, whose preserved rows ``_reflect_trace`` then discarded.
            # ``ZeroOperator``'s symmetric space hooks are the designed
            # mechanism for a genuine map BETWEEN spaces: the forward emits the
            # zero of Γ₋, the transpose the zero of Γ₊. Relying on the
            # endomorphic ``0.0 * x`` echo would be wrong in principle and
            # merely lucky in practice (|Γ₊| = |Γ₋| on every reachable
            # fixture, so the shapes coincide — an accident, not a contract).
            gamma_out = _outflow_restriction(method_space, "vacuum")
            return stamp_boundary_role(
                _narrowed_zero_operator(
                    method_space, gamma_out, law_key="vacuum",
                )
            )

        if isinstance(law, ReflectiveBoundary):
            # B3.2 — the specular mirror, narrowed to Γ₊ → Γ₋. The pairing
            # itself is built by ``_specular_kernel``, which the albedo arm
            # below also reaches: the mirror sits in this law's G (a symmetry
            # of the domain) and in that law's R (a polished wall), and they
            # realize to the same matrix.
            gamma_out = _outflow_restriction(method_space, "reflective")
            return _attenuated_kernel_operator(
                _specular_kernel(
                    quad, method_space, gamma_out,
                    axis=law.axis, law_key="reflective",
                ),
                law.albedo,
                method_space=method_space, gamma_out=gamma_out, law_key="reflective",
            )

        if isinstance(law, WhiteBoundary):
            # B3.4a — the Lambertian kernel, narrowed to Γ₊ → Γ₋.
            #
            # the Lambertian arm derives BOTH half-traces from the one
            # face-name → signed-projection primitive, classified against
            # ``TANGENTIAL_EPS`` — so its old private ``> 0.0`` outflow test is
            # gone, and with it the disagreement with the trace space on every
            # quadrature carrying tangential ordinates.
            gamma_out = _outflow_restriction(method_space, "white")
            return _attenuated_kernel_operator(
                _checked_angular_average(
                    quad, method_space, gamma_out,
                    axis=law.axis, outward_sign=law.outward_sign,
                    law_key="white",
                ),
                law.albedo,
                method_space=method_space, gamma_out=gamma_out, law_key="white",
            )

        if isinstance(law, AlbedoBoundary):
            # B3.4b — a SURFACE returning α of the flux, in the angular shape
            # its closure names. ``G`` is ``SelfPairedDeck.identity()``
            # unconditionally (a wall is not a quotient of the domain), so
            # the whole law is in
            # ``R`` and both arms land in the same two bodies the
            # geometry-tier laws use.
            #
            # Dispatch reads the FACTOR, not ``law.reemission``. The affine
            # form's tier is the public surface the campaign built; the
            # closure is the field the law happens to store, and reaching for
            # it would bypass the abstraction. Reading ``response_kernel``
            # also single-sources α (``kernel.amplitude``) and is the shape
            # B4 generalizes: that phase turns this branch's isinstance chain
            # into the realizer's ONE dispatch, over all laws.
            gamma_out = _outflow_restriction(method_space, "albedo")
            kernel = law.response_kernel
            if isinstance(kernel, SpecularReemission):
                return _attenuated_kernel_operator(
                    _specular_kernel(
                        quad, method_space, gamma_out,
                        axis=kernel.axis, law_key="albedo",
                    ),
                    kernel.amplitude,
                    method_space=method_space, gamma_out=gamma_out,
                    law_key="albedo",
                )
            if isinstance(kernel, LambertianReemission):
                return _attenuated_kernel_operator(
                    _checked_angular_average(
                        quad, method_space, gamma_out,
                        axis=kernel.axis, outward_sign=kernel.outward_sign,
                        law_key="albedo",
                    ),
                    kernel.amplitude,
                    method_space=method_space, gamma_out=gamma_out,
                    law_key="albedo",
                )
            if isinstance(kernel, ScalarResponse) and kernel.is_zero:
                # α = 0 is NOT under-determined — nothing returns, so no
                # pairing is needed and every closure would agree. It is
                # refused for the other reason: it is a VacuumInflow twin, and
                # admitting it would give one physical law two spellings with
                # two realization paths. Refusing also keeps SN's albedo
                # admission a single uniform rule ("state a closure"), rather
                # than one that turns on an exact float compare — a law that
                # realizes at α = 0 and refuses at α = 1e-300 is a worse
                # contract than one that refuses both.
                raise BoundaryError(
                    f"SNBoundaryRealizer will not realize AlbedoBoundary"
                    f"(albedo=0.0) with no re-emission closure: a surface "
                    f"that returns nothing IS a vacuum, and VacuumInflow() "
                    f"already spells it — with a narrowed zero map that "
                    f"carries both space hooks and an honest transpose. Use "
                    f"VacuumInflow(), or state a closure if you mean a "
                    f"perfectly absorbing wall that should still certify its "
                    f"pairing (AlbedoBoundary(0.0, SpecularReturn(axis=...)) "
                    f"realizes, and to the same zero map).",
                    law="albedo",
                )
            if isinstance(kernel, ScalarResponse):
                raise BoundaryError(
                    f"SNBoundaryRealizer cannot realize AlbedoBoundary"
                    f"(albedo={law.albedo}) with no re-emission closure: on "
                    f"an ANGULAR trace the law is under-determined. Its "
                    f"response R = α·I is an endomorphism of Γ₊ and its "
                    f"geometry G = SelfPairedDeck.identity() supplies no crossing, "
                    f"so "
                    f"nothing says which outgoing direction feeds which "
                    f"incoming one — composing it anyway would pair them by "
                    f"ARRAY POSITION, an artefact of index order carrying no "
                    f"geometry. State the shape: "
                    f"AlbedoBoundary({law.albedo}, SpecularReturn(axis=...)) "
                    f"for a polished surface, or "
                    f"AlbedoBoundary({law.albedo}, IsotropicReturn(axis=..., "
                    f"outward_sign=...)) for a diffuse one. (A SCALAR method "
                    f"needs no closure — α alone is the complete partial-"
                    f"current law, which is why the diffusion realizer takes "
                    f"this same object unchanged.)",
                    law="albedo",
                )
            raise BoundaryError(
                f"SNBoundaryRealizer cannot realize AlbedoBoundary whose "
                f"response kernel is {kernel!r}: no realization body for that "
                f"kernel. Realizable kernels: SpecularReemission, "
                f"LambertianReemission.",
                law="albedo",
            )

        if isinstance(law, PeriodicBoundary):
            # B3.4c — the torus quotient, narrowed to Γ₊(partner) → Γ₋(face).
            #
            # Periodic is the only law whose DOMAIN is not the face it is
            # installed on: ``γ₋ψ|_f = γ₊ψ|_{f'}``. Which face f' is depends on
            # where the law sits, so the geometry factor is asked
            # (``G.domain_face``) rather than the partner being stored — and
            # that call is also where a wrap declared for the wrong axis is
            # refused. The composition supplies the partner's half-trace;
            # ``SNBoundaryOperator._face_domains`` is the consumer.
            _assert_wrap_identification(law.geometry_map, method_space)
            # The body is the IDENTITY on the local index, and that is EARNED,
            # not assumed: the guard above proves Γ₊(partner) and Γ₋(face) are
            # the same global ordinate set, so ordinate n of the partner's
            # outflow lands at ordinate n of this face's inflow with no
            # relabelling. (It is a geometric theorem — the outward normals are
            # opposite, so ``n_f = −n_f'`` sends outgoing to incoming — but a
            # theorem the realization checks rather than trusts.)
            return stamp_boundary_role(IdentityOperator() & IdentityOperator())

        if isinstance(law, PrescribedInflow):
            # P3 — realize the LINEAR factor, which is ZERO. The law is
            # affine, ``γ₋ψ = L γ₊ψ + q``, and this tier realizes ``L``
            # alone; the source ``q`` travels the boundary-source channel
            # (:ref:`bc-affine-source-channel`), assembled from the declared
            # law by ``AngularBoundarySourceSink.from_mesh_laws``. So
            # prescribed inflow returns the SAME expression vacuum does,
            # which is the honest statement of the algebra: the two laws
            # differ only in a term that is not an operator.
            #
            # ⛔ Until P3 this arm returned an ``IncomingSourceOperator``
            # whose ``apply`` IGNORED its input and emitted ``q`` — an
            # AFFINE map in a linear slot. It carried no BOUNDARY stamp, and
            # the stamp's absence was believed to fence it out of ``B``. It
            # did not: ``SNBoundaryOperator._face_laws`` collects EVERY
            # face's law with no role filter, so the source was delivered
            # through ``B`` as well. MEASURED on a declared inflow at
            # ``q = 2.5``: ``|B(0)| = 2.5`` and ``|B(2x) − 2B(x)| = 2.5``,
            # and once P2′ wired the source channel the converged flux
            # carried the inflow TWICE (ratio 2.000000 against the vacuum
            # control). Returning the zero morphism is what makes the two
            # channels one.
            #
            # Prescribed inflow ignores its input, but a law that cannot
            # name its own domain is not realized, it is guessed — so the
            # Γ₊ restriction is still built, now for its VALUE as well as
            # its guard: the zero map is bound with it.
            gamma_out = _outflow_restriction(method_space, "prescribed_inflow")
            return stamp_boundary_role(
                _narrowed_zero_operator(
                    method_space, gamma_out, law_key="prescribed_inflow",
                )
            )

        raise BoundaryError(
            f"SNBoundaryRealizer cannot realize "
            f"{type(law).__name__} — no isinstance dispatch case. "
            f"Available cases: VacuumInflow, "
            f"ReflectiveBoundary, WhiteBoundary, "
            f"AlbedoBoundary, PeriodicBoundary, "
            f"PrescribedInflow. For rank-N compositions built via "
            f"Wave-0 OperatorSum/ScaledOperator algebra over "
            f"BoundaryTraceLaw leaves, use "
            f"orpheus.geometry.boundary.realize_recursively with "
            f"realizer=SNBoundaryRealizer().",
            law=type(law).__name__,
        )


# The rank-N composition walker (``realize_recursively``) lived here
# until #290 P7b. With the second functional realizer shipped
# (diffusion, #290 P3) the walk's realizer-genericity became real and
# the walker moved to its method-blind home,
# :func:`orpheus.geometry.boundary.realize_recursively` — the leaf
# realizer is now a REQUIRED argument (``SNBoundaryRealizer()`` for the
# SN path). Single BCs never routed through it: ``SNMesh``'s
# ``realize_boundary_law`` arm (the ``TransportMethod`` hook) calls
# ``SNBoundaryRealizer().realize`` directly.

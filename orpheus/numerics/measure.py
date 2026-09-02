r"""Discrete measures and quadrature composition primitives.

A *quadrature rule* is mathematically a **discrete measure** on a
measurable space :math:`(\mathcal{X}, \Sigma)` — a finite weighted
sum of Dirac atoms,

.. math::
   :label: discrete-measure-definition

   \mu = \sum_{i=1}^{N} w_i \, \delta_{x_i},
   \qquad x_i \in \mathcal{X}, \quad w_i \in \mathbb{R}.

The Lebesgue integral of an integrable test function :math:`f` against
:math:`\mu` is the familiar quadrature sum,

.. math::
   :label: discrete-measure-integrate

   \int_{\mathcal{X}} f \, d\mu \;=\; \sum_{i=1}^{N} w_i \, f(x_i).

This module promotes that view to a first-class primitive. The
project's single quadrature class
(:class:`~orpheus.numerics.quadrature.Quadrature`, reached through its
named factories) composes 1-D rules into 2-D / :math:`S^2` rules
**internally**, hiding the tensor-product structure behind
:meth:`~orpheus.numerics.quadrature.Quadrature.product`. With
:class:`DiscreteMeasure` the four canonical operations are exposed:

- **Tensor product** ``μ * ν``. Product measure on
  :math:`\mathcal{X} \times \mathcal{Y}` — Fubini-Tonelli on discrete
  factors, in place of a product implicit inside a rule constructor.
- **Direct sum** ``μ + ν``. Concatenation on a shared space — used
  when a domain is partitioned into subintervals each carrying its
  own rule.
- **Pushforward** ``μ.pushforward(φ)``. Image measure
  :math:`\varphi_* \mu` under a measurable map
  :math:`\varphi : \mathcal{X} \to \mathcal{Y}`. The change-of-variables
  identity holds verbatim:
  :math:`\int g \, d(\varphi_* \mu) = \int g \circ \varphi \, d\mu`.
- **Restriction** ``μ.restrict(E)``. Indicator-multiplication
  :math:`\mathbf{1}_E \cdot \mu` — used for half-range SN sweeps and
  for cutting bundles to a region.

The :class:`BundleMeasure` class implements **disintegration**: a
measure on a fibered space :math:`\pi : \mathcal{X} \to \mathcal{B}`
is given by a base measure on :math:`\mathcal{B}` plus, for each base
point, a fiber measure on :math:`\pi^{-1}(b)`. This is the structure
MoC ray quadratures need (parallel ray bundles per polar angle), and
is shipped here so the abstraction is correct from day one — the SN
side of the campaign does not consume :class:`BundleMeasure` directly,
but the MoC migration in Wave 2 does.

References
----------

* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. Springer. Chapter 3 (numerical integration);
  Theorem 3.6.20 (Gauss-Legendre exactness for polynomials of
  degree :math:`\le 2N-1`).
* Trefethen, L.N. (2008). "Is Gauss quadrature better than
  Clenshaw-Curtis?" *SIAM Review* **50**, 67-87. Measure-theoretic
  perspective on spectral integration rules.
* Bourbaki, N. (1969). *Éléments de mathématique. Intégration*,
  Chapter VI §3 (disintegration of measures).

See Also
--------

:ref:`discrete-measures` (theory page) — composition algebra,
metadata-propagation table, structural connections to
:class:`~orpheus.numerics.quadrature.Quadrature` and the upcoming
``invariance_group`` tag from Issue 3.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Callable, Literal, Protocol, runtime_checkable

import numpy as np

# ⭐ A module-scope import, and it is safe by design rather than by luck:
# ``manifold`` reaches ``symmetry`` (which imports THIS module at module scope)
# only under ``TYPE_CHECKING``, so the cycle ``measure → manifold → symmetry →
# measure`` never closes at runtime. ``[M]`` 2026-09-01, injected and run in
# five import orders (measure / symmetry / geometry / registry / sn.solver
# first) — all clean. Do not promote ``manifold``'s TYPE_CHECKING import of
# ``symmetry`` to module scope without re-running that probe.
from orpheus.numerics.manifold import (
    COSINE_INTERVAL,
    EnergyGroups,
    Interval,
    Manifold,
    Quotient,
    RealSpace,
    Sphere,
)

if TYPE_CHECKING:
    # ``axis`` is lower-level (imports nothing from this module); the
    # forward-ref keeps the annotation cost-free — the runtime import
    # lives inside the ``axis()`` mint, mirroring ``space``'s lazy style.
    from orpheus.numerics.axis import Axis

    # Forward-ref only: ``symmetry`` imports FROM this module, so a runtime
    # import would cycle. The field stores a ``SubgroupOfO3`` object (passed
    # by the angular quadrature factories, which already import it).
    from orpheus.numerics.symmetry import SubgroupOfO3
    from orpheus.numerics.space import FunctionSpace
    # Same forward-ref reason: ``generating_measure`` imports FROM this
    # module (its ``gauss`` returns a DiscreteMeasure), so the arrow runs
    # one way and the field holds the object by forward reference.
    from orpheus.numerics.generating_measure import GeneratingMeasure
    # ``exactness`` imports ``DiscreteMeasure`` from here, so the same
    # one-way arrow applies: the field holds an ``ExactnessClaim`` by
    # forward reference.
    from orpheus.numerics.exactness import ExactnessClaim

# ---------------------------------------------------------------------------
# The point set a measure lives on
# ---------------------------------------------------------------------------

# ⛔ Until 2026-09-01 this section defined ``Space = str`` plus six ``SPACE_*``
# string tags, and its comment read: *"Treated as opaque labels at runtime —
# composition operations stitch them together (e.g. `"[-1,1]" * "S^1"` →
# `"[-1,1] × S^1"`) but no type-level enforcement is attempted, in line with
# the design note in `.claude/plans/sn_reshape.md` Issue 2 ('Don't try to
# enforce ``Space`` types via Python generics — not expressive enough without
# runtime overhead')."*
#
# ⭐ That note was answered rather than overruled. The objection was to
# *generics* — a phantom type parameter, erased at runtime, buying nothing;
# it is still correct. What ships instead is a closed sum of ORDINARY VALUES
# (:mod:`orpheus.numerics.manifold`), so the "stitching" the comment describes
# is now :meth:`Manifold.__mul__` returning a ``Product`` that can answer
# ``dim`` and ``contains``, at no runtime cost and with no generics.
#
# ``Space = str`` was retired 2026-09-01 (#429 tracker 2.0c): a measure now
# names the point set it lives on with a :class:`~orpheus.numerics.manifold.Manifold`,
# so ``support`` carries the object's structure — its dimension, its membership
# predicate, its quotient and product algebra — instead of a tag that only a
# reader could interpret. The six ``SPACE_*`` string constants below retired
# with it; their replacements are the shipped ``Manifold`` members
# (``REAL_LINE``, ``HALF_LINE``, ``COSINE_INTERVAL``, ``UNIT_INTERVAL``,
# ``CIRCLE``, ``SPHERE``, ``ENERGY``) and every one of them reproduces the
# retired tag's spelling bit-identically through :attr:`Manifold.name`
# (``[M]`` 10 of 10 at the retype).
# ⭐ The circle names the MANIFOLD, not a chart of it — a distinction that
# outlived the tag it was written for, and is now carried by the type. The tag
# was ``"[0,2π)"`` until 2026-08-02, which made it assert a COORDINATE; the
# load-bearing consequence is that the periodic trapezoid's nodes are the roots
# of unity as POINTS (:func:`~orpheus.numerics.roots_of_unity.roots_of_unity`),
# because only in that representation is an on-axis node exactly ``0.0`` — an
# angle chart would have re-introduced the very rounding
# :mod:`~orpheus.numerics.roots_of_unity` exists to remove. Uniform measure on
# S^1 has mass 2π, on S^2 mass 4π — one family, one rule. That is now
# :class:`~orpheus.numerics.manifold.Circle` and
# :class:`~orpheus.numerics.manifold.Sphere`, which cannot be spelled as a chart.

# The physical FACTOR of transport phase space (position × direction × energy)
# a measure discretises. A closed enumeration — there are exactly these factors
# (steady-state transport) — so :attr:`DiscreteMeasure.phase` is a derived
# *category*, NOT an arbitrary string. The three factors are genuinely
# different objects (a compact O(3)-sphere vs an unbounded energy half-line vs
# a bounded mesh domain); each earns its own typed support as it is first
# Frame-bound (the angular factor is the worked instance — see :attr:`phase`).
Phase = Literal["angular", "spatial", "energy"]


# ---------------------------------------------------------------------------
# Core primitive
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DiscreteMeasure:
    r"""A finite weighted point measure on a measurable space.

    Algebraic shape: :math:`\mu = \sum_i w_i \, \delta_{x_i}`.
    See :eq:`discrete-measure-definition`.

    Parameters
    ----------
    nodes : np.ndarray
        Array of node coordinates. Shape is one of:

        - ``(N,)`` — 1-D measure (a 1-D node is a scalar).
        - ``(N, d)`` — d-dimensional measure with ``N`` atoms in
          :math:`\mathbb{R}^d`. ``d`` is the geometric dimension of
          ``space``; tensor-product measures stack dimensions.

        Conversion from ``(N,)`` to ``(N, 1)`` is at the user's
        discretion — both shapes are accepted and routed correctly
        in :meth:`integrate`. Composition operations preserve the
        caller's choice.
    weights : np.ndarray
        Array of weights, shape ``(N,)``. Weights MAY be negative
        (e.g. signed quadratures) but most rules in this codebase
        produce non-negative weights.
    support : Manifold
        The **point set the atoms live on** — :class:`Sphere`,
        :class:`Interval`, :class:`RealSpace`, an orbit space
        :class:`Quotient`, … (:mod:`orpheus.numerics.manifold`).

        It is an OBJECT, not a tag, and three things are derived from it that
        a string could only assert: :attr:`phase` (which factor of transport
        phase space this discretises) reads its TYPE; :attr:`space`'s name
        reads its :attr:`~Manifold.name`; and :meth:`quotient` folds onto
        ``support.quotient(G)`` — the catalogue's orbit space, not a
        fabricated one. Direct sum requires equal supports and the tensor
        product multiplies them, both through the manifold's own algebra.

        NOTE: this is the *continuous support*, distinct from the induced
        *discrete function space* :attr:`space` (the operator-algebra domain)
        and from the physical :attr:`phase` (angular/spatial/energy).
    invariance_group : SubgroupOfO3 | None, optional
        The maximal subgroup of :math:`O(3)` under which the measure is
        invariant (its nodes permute among themselves with weights
        preserved), as a typed
        :class:`~orpheus.numerics.symmetry.SubgroupOfO3` — e.g.
        ``SubgroupOfO3.OctahedralOh`` for a Lebedev cubature,
        ``SubgroupOfO3.Dnh(n_phi)`` for a product rule, ``SubgroupOfO3.Mirror('x')``
        for the slab's polar marginal — never a continuous group for a finite
        point set, which is closed under :math:`SO(2)_a` only on the axis
        (ERR-072; the slab's rule *spent* ``SO2('x')`` instead, see
        :attr:`quotient_group`). ``None`` means "unspecified". This symmetry is the structural
        signature of the **angular** phase-space factor: a measure carrying
        an :math:`O(3)`-subgroup invariance lives on :math:`S^2` (Erlangen —
        the group fixes the homogeneous space), which :attr:`phase` reads.
    exactness : ExactnessClaim | None, optional
        What this rule integrates exactly, stated **whole**: the
        reference measure the claim is against, and the degree indexing
        that reference's own orthogonal system
        (:class:`~orpheus.numerics.exactness.ExactnessClaim`).

        .. math::

           \sum_i w_i \, f(x_i) \;=\; \int_\mathcal{X} f \,
           d\lambda \qquad \text{for every } f \in V.

        ``None`` if the rule has no exactness claim or it has not been
        computed.

        **This used to be two loose fields** — a bare
        ``degree_of_exactness`` beside a ``generating_measure`` — and a
        degree could therefore exist without saying what it was about.
        Two measured consequences, and they are *different* bugs:

        * same degree, different **measure**: `[M]` Gauss-Chebyshev is
          exact to :math:`2n-1` against
          :math:`\int q\,(1-x^2)^{-1/2}dx` and misses
          :math:`\int q\,dx` by **0.696** at :math:`n=4,\ q=x^6`;
        * same degree, different **space**: `[M]` :math:`n` equispaced
          nodes on an interval are the midpoint rule (*algebraic*
          degree 1); the same nodes on the circle are the periodic
          trapezoid (*trigonometric* degree :math:`n-1`).

        :attr:`degree_of_exactness` and :attr:`generating_measure`
        survive as **derived read-only views** over this one field, so
        the convenient readings remain available while the storage that
        let them drift apart is gone.

        A claim whose reference is a
        :class:`~orpheus.numerics.generating_measure.GeneratingMeasure`
        **cannot over-claim** — the degree follows from the
        construction. A rule without one rests on external authority (a
        published table plus a citation, as Lebedev and level-symmetric
        :math:`S_N` do). `[M]` That is exactly how issue #327 was
        possible: ``level_symmetric_sn`` hand-assigned one equal weight
        to every ordinate, so nothing constrained its advertised
        degree, and it WAS degree-3 at every order while claiming
        :math:`N-1` — until 2026-08-06, when the weights became
        per-:math:`O_h`-orbit moment-solved and the gate
        (``tests/numerics/test_advertised_degree_is_measured.py``)
        began asserting measured == advertised both directions.

    Notes
    -----

    The class is **frozen** (immutable). Composition operations return
    new instances; mutation is not supported. This matches the
    mathematical convention that a measure is a fixed object, not a
    state to be evolved.

    Composition algebra (see :ref:`discrete-measures` for the full
    propagation table):

    - ``μ * ν`` — tensor product, returns a measure on
      ``f"{μ.support} × {ν.support}"``.
    - ``μ + ν`` — direct sum on a shared space; raises ``ValueError``
      if ``μ.support != ν.support``.
    - ``μ.pushforward(φ)`` — image measure under ``φ``. Weights are
      preserved; nodes become ``φ(nodes)``. The Jacobian of ``φ`` is
      the **caller's responsibility**: this is the φ-image semantics,
      not a Radon-Nikodym derivative against a reference measure.
    - ``μ.quotient(G)`` — the fold: one atom per :math:`G`-orbit on
      ``f"{μ.support}/{G.name}"``, carrying the summed orbit weight
      :math:`W = w \cdot |G| / |\mathrm{Stab}|`. Mass is preserved;
      refuses unless the measure is certified :math:`G`-invariant.
    - ``μ.restrict(E)`` — keep atoms where ``E(x)`` is true; drop the
      rest. Weights of kept atoms are preserved.
    """

    nodes: np.ndarray
    weights: np.ndarray
    support: Manifold
    invariance_group: "SubgroupOfO3 | None" = None
    exactness: "ExactnessClaim | None" = None

    def __post_init__(self) -> None:  # pragma: no cover - guard
        # Invariants enforced at construction. Frozen dataclasses
        # forbid attribute mutation, so we never re-check.
        if self.weights.ndim != 1:
            raise ValueError(
                f"weights must be 1-D, got shape {self.weights.shape}"
            )
        if self.nodes.shape[0] != self.weights.shape[0]:
            raise ValueError(
                f"nodes and weights disagree on N: nodes.shape="
                f"{self.nodes.shape}, weights.shape={self.weights.shape}"
            )
        if self.nodes.ndim not in (1, 2):
            raise ValueError(
                f"nodes must be 1-D (N,) or 2-D (N, d), got shape "
                f"{self.nodes.shape}"
            )

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def n_points(self) -> int:
        """Number of atoms (length of ``weights``)."""
        return int(self.weights.shape[0])

    @property
    def dim(self) -> int:
        """Geometric dimension :math:`d` of the node array.

        Returns ``1`` for ``(N,)``-shaped nodes, otherwise ``nodes.shape[1]``.
        """
        if self.nodes.ndim == 1:
            return 1
        return int(self.nodes.shape[1])

    @property
    def space(self) -> "FunctionSpace":
        r"""The induced discrete-:math:`L^2` :class:`FunctionSpace` — the domain.

        A discrete measure :math:`\mu = \sum_n w_n \delta_{x_n}` induces an
        :math:`L^2(\mu)` space of sampled values: per-node arrays of shape
        ``(N,)`` under the inner product
        :math:`\langle f, g\rangle_\mu = \sum_n w_n f_n g_n` (the weights ARE
        the metric). This is the operator-algebra **domain** a
        :class:`~orpheus.numerics.frame.FrameBase` reads — mirroring the way a
        :class:`~orpheus.numerics.basis.Basis` owns its coefficient
        :attr:`~orpheus.numerics.basis.Basis.space` (the codomain): the measure
        owns its sampled-value space, so the ``(domain, codomain)`` pair is
        symmetric and neither is fabricated ad-hoc by the binder.

        Distinct from :attr:`support` (the *continuous* measurable space the
        nodes live on, a :class:`~orpheus.numerics.manifold.Manifold`) — this is
        the *discrete* function space.
        Lazy import (``FunctionSpace`` is a higher-level construct; mirrors
        :attr:`SphericalHarmonicBasis.space`'s lazy import).
        """
        from orpheus.numerics.space import FunctionSpace

        return FunctionSpace(
            name=f"L2[{self.support.name}]",
            shape=(self.n_points,),
            inner_product_weights=self.weights,
        )

    def axis(self, label: str) -> "Axis":
        r"""Mint the space-factor :class:`~orpheus.numerics.axis.Axis` this
        measure generates — the axis-composed sibling of :attr:`space`.

        An axis is a **forgetful map** from its generator: it keeps the
        weights (the factor measure) and drops the nodes. Minting through
        the generator makes the forgetting recoverable — the returned
        axis carries ``generator=self`` (provenance, excluded from axis
        identity), so a consumer can read the un-forgotten data through
        the space instead of being handed the measure separately (CS5).

        ``kind`` is NOT a parameter: a discrete measure's samples are
        point values with a coordinate cone, so the generator's type
        implies :attr:`~orpheus.numerics.axis.BasisKind.NODAL` — the
        "nodal basis carrying no nodes" mis-declaration is unspellable
        on this path. (A MODAL axis is minted by its
        :class:`~orpheus.numerics.basis.Basis`, not by a measure.)

        Lazy import for style symmetry with :attr:`space`; the runtime
        edge ``measure → axis`` is new and acyclic (``axis`` imports
        nothing from this module — [M] 2026-08-29, CS5 ground memo §F).
        """
        from orpheus.numerics.axis import Axis, BasisKind

        return Axis(
            label,
            (self.n_points,),
            weights=self.weights,
            kind=BasisKind.NODAL,
            generator=self,
        )

    @property
    def quotient_group(self) -> "SubgroupOfO3 | None":
        r"""The group this measure's support was folded BY, or ``None``.

        ⭐ **Derived, never stored** — tracker 2.0d asked for a
        ``quotient_group`` FIELD and it dissolved into 2.0c at the phase
        opener: :attr:`Quotient.by` already *is* the group, so a stored copy
        would be a second home for one fact and the two would have to be kept
        in agreement by hand across every operation that rebuilds a measure
        (:meth:`restrict`, :meth:`consolidate`, :meth:`partition_by`,
        :meth:`on_orbit_space`, …).

        ⚠ Distinct from :attr:`invariance_group`, and the two are almost
        complementary: ``invariance_group`` is the symmetry the measure
        **HAS**, ``quotient_group`` the symmetry it has **SPENT**. A folded
        rule reports ``None`` for the first (the mirror was quotiented away,
        so its nodes are no longer closed under it) and the group for the
        second. Reading either as the other is
        ``plan-authoring`` §3's ambiguous-name hazard, and ``invariance_group``
        has drawn it before.
        """
        return self.support.by if isinstance(self.support, Quotient) else None

    @property
    def phase(self) -> Phase:
        r"""Which factor of transport phase space this measure discretises.

        Transport lives on phase space = position × direction × energy. A
        measure always discretises exactly one factor; :attr:`phase` is *which*
        — a derived **category**, read from the measure's structure rather than
        a hand-set label (see :data:`Phase`).

        Derivation, by category — a ``match`` on the TYPE of the support
        manifold (the user's ruling of 2026-09-01: ``Phase`` is a transport
        taxonomy and a manifold is pure geometry, so the dispatch lives on
        the measure and reads the manifold's type):

        * **angular** — iff :attr:`support` is a :class:`Sphere` or any
          quotient of one: Lebedev's :math:`S^2`, the fold's
          :math:`S^2/\sigma_y`, the slab's :math:`S^2/SO(2)_x` (declared
          by :func:`~orpheus.numerics.quadrature.rules_1d.gauss_legendre_on_polar_orbit`
          since tracker 2.4). The direction variable lives on the sphere,
          and every residual symmetry a rule spent still leaves it there.
          ⚠ One more input is angular by a *fallback*: a rule on a bare
          :class:`Interval` that carries an :math:`O(3)` invariance tag —
          the chart-level :math:`\mu`-rule
          :func:`~orpheus.numerics.quadrature.rules_1d.gauss_legendre_on_mu`
          itself, which serves as the polar FACTOR of product rules and so
          cannot declare an orbit space. It is the one shipped input that
          reaches that arm (`[M]` 2026-09-01).
        * **spatial** — iff :attr:`support` is a :class:`RealSpace`. Spatial
          measures carry no :math:`O(3)` symmetry (a mesh is not
          rotation-invariant); they are a *different category*.
        * **energy** — iff :attr:`support` is :class:`EnergyGroups`. The
          multigroup energy axis: a **counting** measure on the discrete group
          index (``weights = 1``; ``φ_g`` is group-integrated), Frame-bound by
          energy condensation
          (:meth:`~orpheus.data.macro_xs.mixture.Mixture.condense`, via
          :meth:`~orpheus.data.energy_grid.EnergyGrid.as_measure`).

        Anything else raises — an untagged rule on an ambiguous point set
        has no determined phase. The asymmetry IS the signal: each factor
        gets its own typed machinery as it is bound, not a premature uniform
        abstraction (the user's design ruling). The phase cannot be read off
        the bare nodes — a :math:`\mu\in[-1,1]` chart is geometrically
        indistinguishable from a 1-D spatial interval, and an energy group
        index from any integer-noded counting rule; the physical identity is
        exactly what the support manifold supplies.
        """
        match self.support:
            # A sphere and any quotient of one are the direction variable's
            # home, whatever residual symmetry the rule spent. The quotient arm
            # is what the shipped fold needs: ``folded_product(4, 8)`` lives on
            # S²/σ_y with ``invariance_group=None`` (the mirror was quotiented
            # AWAY, so the measure is not invariant under it), and before this
            # arm existed asking it for its phase RAISED.
            case Sphere() | Quotient(base=Sphere()):
                return "angular"
            case EnergyGroups():
                return "energy"
            case RealSpace():
                return "spatial"
            case _ if self.invariance_group is not None:
                # ⚠ The one genuinely ambiguous point set, and the reason
                # ``phase`` is a property of the MEASURE and not of the
                # manifold: a rule on ``Interval(-1, 1)`` spells itself
                # exactly as a 1-D spatial interval would. The O(3) symmetry
                # is what tells them apart, and it lives on the measure.
                #
                # Since tracker 2.4 (2026-09-01) the SLAB no longer needs
                # this arm — its quadrature declares S^2/SO2_x and answers
                # from the Sphere-quotient arm above. ⛔ The arm is NOT
                # unreachable, though, which the 2.4 pre-flight predicted
                # it would be: the chart-level μ-rule `gauss_legendre_on_mu`
                # keeps a bare Interval with a Mirror('x') tag, because it is
                # also the polar FACTOR of every product rule and a factor
                # declaring an orbit space about x inside a product about z
                # would be a false claim. `[M]` that rule is the one shipped
                # input still arriving here.
                return "angular"
        raise NotImplementedError(
            f"phase is undetermined for a measure on {self.support.name!r}: "
            f"the direction factor is any sphere or quotient of one, the "
            f"energy factor is EnergyGroups, and the position factor is "
            f"RealSpace. An Interval is ambiguous between the slab's μ-axis "
            f"and a 1-D spatial interval, so it is angular only when the "
            f"measure carries the O(3) symmetry that says so."
        )

    # ------------------------------------------------------------------
    # Lebesgue integral against the discrete measure
    # ------------------------------------------------------------------

    def integrate(
        self,
        f: Callable[[np.ndarray], np.ndarray] | np.ndarray,
    ) -> np.ndarray:
        r"""Compute :math:`\int f \, d\mu \;=\; \sum_i w_i \, f(x_i)`.

        See :eq:`discrete-measure-integrate`.

        Parameters
        ----------
        f : callable | np.ndarray
            Either a vectorised test function called once with the
            full ``nodes`` array (returning an array of shape ``(N,)``
            or ``(N, k)`` for vector-valued integrands), OR a
            pre-evaluated value array of shape ``(N,)`` or ``(N, k)``.
            The array overload is the load-bearing affordance for
            ``mesh.volume_measure(scalar_flux)`` — at production
            integration sites where ``f(nodes)`` would be redundant
            (the values are already in hand), passing the array
            directly avoids the round-trip through a ``lambda x:
            values`` wrapper.

        Returns
        -------
        np.ndarray
            ``Σ_i w_i f(x_i)``. Shape is ``()`` for scalar
            integrands or ``(k,)`` for vector-valued integrands.

        Notes
        -----

        When ``f`` is callable, it is called **once** with the full
        node array, not per-node. Vectorise your integrand (use
        ``np.cos``, not ``math.cos``) — the project's test suite
        assumes this.
        """
        if callable(f):
            values = f(self.nodes)
        else:
            values = f
        values_arr = np.asarray(values)
        if values_arr.ndim == 1:
            return np.dot(self.weights, values_arr)
        # Vector-valued integrand: contract the leading axis.
        return np.einsum("i,i...->...", self.weights, values_arr)

    # ------------------------------------------------------------------
    # Dunder ergonomics — measure as functional, iterable, indexed
    # ------------------------------------------------------------------

    def __call__(
        self,
        f: Callable[[np.ndarray], np.ndarray] | np.ndarray,
    ) -> np.ndarray:
        """Alias for :meth:`integrate`. Lets ``mu(f)`` read as math."""
        return self.integrate(f)

    def __iter__(self):
        """Iterate over ``(node, weight)`` pairs."""
        for i in range(self.n_points):
            yield (self.nodes[i], self.weights[i])

    def __len__(self) -> int:
        """Return :attr:`n_points`."""
        return self.n_points

    def __getitem__(self, i: int) -> tuple[np.ndarray, float]:
        """Return the ``i``-th ``(node, weight)`` pair."""
        return (self.nodes[i], self.weights[i])

    def __repr__(self) -> str:
        if self.invariance_group is None:
            return (
                f"DiscreteMeasure(support={self.support.name!r}, "
                f"n_points={self.n_points})"
            )
        return (
            f"DiscreteMeasure(support={self.support.name!r}, "
            f"n_points={self.n_points}, "
            f"invariance_group={self.invariance_group!r})"
        )

    # ------------------------------------------------------------------
    # Tensor product (Fubini-Tonelli on discrete factors)
    # ------------------------------------------------------------------

    def __mul__(self, other: DiscreteMeasure) -> DiscreteMeasure:
        r"""Tensor product :math:`\mu \otimes \nu`.

        Returns the product measure on
        :math:`\mathcal{X} \times \mathcal{Y}`. With :math:`N_1 = |\mu|`
        and :math:`N_2 = |\nu|`, the result has :math:`N_1 N_2` atoms.

        Node convention:

        - If both factors are 1-D (``dim == 1``), the result has
          shape ``(N_1 N_2, 2)``, with column 0 carrying ``μ``-nodes
          and column 1 carrying ``ν``-nodes (Cartesian-product order
          — outer loop ``μ``, inner loop ``ν``, matching
          ``np.meshgrid(..., indexing='ij')``).
        - If either factor is multi-dimensional, both factors are
          first promoted to ``(N, d)`` form via ``[:, None]`` and
          then horizontally stacked; result shape ``(N_1 N_2, d_1
          + d_2)``.

        Weights are the outer product flattened in Cartesian order:
        ``w_ij = μ.weights[i] * ν.weights[j]``.

        Metadata propagation (see :ref:`discrete-measures` for the
        full table):

        - ``space`` becomes ``f"{μ.support} × {ν.support}"``.
        - ``invariance_group`` becomes ``None`` unless both factors
          carry compatible groups (Issue 3 will refine this).
        - ``degree_of_exactness`` becomes ``min(p_μ, p_ν)`` if both are
          set **and the two factors share a generating measure**, else
          ``None``. (For separable polynomials of degree :math:`p` per
          axis, both rules need to integrate monomials up to :math:`p`
          exactly.)
        - ``generating_measure`` becomes ``None``: the product rule's
          reference is the *product* measure :math:`w_\mu \otimes
          w_\nu`, which a 1-D
          :class:`~orpheus.numerics.generating_measure.GeneratingMeasure`
          cannot name.

        .. note:: Why the degree needs the shared-measure condition

           A degree of exactness is a claim relative to a measure, so
           combining two rules with *different* references produces a
           claim about a mongrel measure no consumer expects. `[M]`
           ``gauss_legendre(4) * gauss_chebyshev(4)`` used to advertise
           ``degree_of_exactness = 7`` while integrating the constant
           ``1`` to **6.2832** instead of 4 — Gauss-Chebyshev is not
           exact for the unweighted integral even at degree **0**
           (3.1416 vs 2). The composite was exact for neither factor's
           measure while claiming both their degrees.

           When both factors share a reference the product measure is
           that measure's product power — a canonical object — and
           ``min`` is the right claim about it. Two untagged rules count
           as sharing (both ``None``), which is the pre-2026-08-02
           behaviour for every rule that does not yet declare a measure.
        """
        # Promote 1-D nodes to column form so we can hstack uniformly.
        a = self.nodes if self.nodes.ndim == 2 else self.nodes[:, None]
        b = other.nodes if other.nodes.ndim == 2 else other.nodes[:, None]
        n1, d1 = a.shape
        n2, d2 = b.shape

        # Cartesian product of nodes: outer loop μ, inner loop ν.
        # Equivalent to np.meshgrid(a_idx, b_idx, indexing='ij'),
        # but vectorised through tile/repeat to avoid creating
        # explicit index grids.
        a_rep = np.repeat(a, n2, axis=0)            # (n1*n2, d1)
        b_tile = np.tile(b, (n1, 1))                # (n1*n2, d2)
        new_nodes = np.hstack([a_rep, b_tile])      # (n1*n2, d1+d2)

        # Outer product of weights, flattened in matching order.
        new_weights = np.outer(self.weights, other.weights).reshape(-1)

        # Metadata propagation. ⭐ The product of the two point sets IS the
        # product manifold — ``Manifold.__mul__`` — not a string that spells one.
        new_space = self.support * other.support

        return DiscreteMeasure(
            nodes=new_nodes,
            weights=new_weights,
            support=new_space,
            invariance_group=None,
            # A TENSOR product, so the reference is the PRODUCT measure —
            # not either factor's. ``tensor_with``, not ``combined_with``:
            # the direct sum below lands on the same space and keeps the
            # shared reference, and one helper serving both is how a
            # product came to inherit a factor's reference and claim
            # exactness against a measure it is not exact against.
            exactness=self._tensor_exactness(other),
        )

    # ------------------------------------------------------------------
    # Direct sum (concatenation on a shared space)
    # ------------------------------------------------------------------

    def __add__(self, other: DiscreteMeasure) -> DiscreteMeasure:
        r"""Direct sum :math:`\mu \oplus \nu` on a shared space.

        Concatenates atoms — useful when a domain is partitioned and
        each piece carries its own quadrature rule (e.g. composite
        Simpson, panelled Gauss-Legendre).

        Requires ``self.support == other.support`` and matching node
        dimensionality (``self.dim == other.dim``). Node arrays are
        concatenated along axis 0; weights are concatenated.

        Metadata propagation:

        - ``space`` is preserved (must already match).
        - ``invariance_group`` is set to ``None`` — concatenation
          generally breaks any invariance the factors had.
        - ``degree_of_exactness`` is set to ``min(p_μ, p_ν)`` if both
          are set **and the two pieces share a generating measure**,
          else ``None`` — the composite rule is at most as exact as its
          weakest piece on the union domain, and only if the two pieces
          are exact against the *same* measure (see :meth:`__mul__` for
          the measurement that motivated the condition).
        - ``generating_measure`` is set to ``None``: the union's
          reference is a measure on the union domain, not either
          piece's.

        Raises
        ------
        ValueError
            If the two measures live on different spaces or have
            different node dimensionality.
        """
        if self.support != other.support:
            raise ValueError(
                f"direct sum requires equal spaces, got "
                f"{self.support.name!r} and {other.support.name!r}"
            )
        if self.dim != other.dim:
            raise ValueError(
                f"direct sum requires equal node dimensions, got "
                f"dim={self.dim} and dim={other.dim}"
            )

        new_nodes = np.concatenate([self.nodes, other.nodes], axis=0)
        new_weights = np.concatenate([self.weights, other.weights])

        return DiscreteMeasure(
            nodes=new_nodes,
            weights=new_weights,
            support=self.support,
            invariance_group=None,
            # A direct sum carries NO claim, and that is a correction.
            # ``__add__`` requires equal supports, so summing two rules
            # for λ yields a rule for **2λ** — its reference is λ₁ + λ₂,
            # not the shared λ. Keeping the shared reference would assert
            # exactness against a measure the sum is not exact against,
            # which is exactly the error the product side just fixed.
            # The predecessor kept ``min(p₁, p₂)`` with the reference
            # dropped — a half-claim, and the state this carve removes.
            # Nothing consumes a direct sum's degree, so the honest
            # answer is no claim rather than a sum-measure type built on
            # speculation.
            exactness=None,
        )

    # ------------------------------------------------------------------
    # Shared combination rule
    # ------------------------------------------------------------------

    def _tensor_exactness(
        self, other: DiscreteMeasure
    ) -> "ExactnessClaim | None":
        r"""The claim of a TENSOR PRODUCT of the two, or ``None``.

        Separate from :meth:`_combined_exactness` because the two
        operations land on **different spaces**, and therefore make
        claims about different reference measures. Until 2026-08-02 one
        helper served both, which is how a product could report a factor's
        reference — asserting exactness against a measure it is not exact
        against. `[M]` a test caught it on the first run of the carve:
        ``gauss_legendre(3) * gauss_legendre(4)`` reported ``legendre``
        where the reference is ``legendre ⊗ legendre``.
        """
        if self.exactness is None or other.exactness is None:
            return None
        return self.exactness.tensor_with(other.exactness)

    # ------------------------------------------------------------------
    # Derived views over the exactness claim
    # ------------------------------------------------------------------

    @property
    def degree_of_exactness(self) -> int | None:
        """The claim's degree, or ``None``.

        A **derived view**, not storage — reading it without
        :attr:`exactness` is reading half a claim, which is precisely
        the state that let a Gauss-Legendre and a Gauss-Chebyshev degree
        look interchangeable. Kept because the degree alone is genuinely
        the useful reading at a call site that already knows the
        reference.
        """
        return None if self.exactness is None else self.exactness.degree

    @property
    def generating_measure(self) -> "GeneratingMeasure | None":
        r"""The claim's reference **when it is a Golub-Welsch generator**.

        A derived view over :attr:`exactness`, and deliberately narrower
        than the claim's own ``reference``: not every reference measure
        generates its rule. The circle's Fourier system has no three-term
        recurrence, and Lebedev's reference (uniform on :math:`S^2`)
        generates nothing at all — both are legitimate references and
        neither is a generating measure, so this returns ``None`` for
        them.

        ``None`` therefore keeps its original meaning: *this rule was not
        built by the Golub-Welsch construction*, and so its exactness
        rests on external authority rather than following from how it was
        made.
        """
        from orpheus.numerics.generating_measure import GeneratingMeasure

        if self.exactness is None:
            return None
        reference = self.exactness.reference
        return reference if isinstance(reference, GeneratingMeasure) else None

    # ------------------------------------------------------------------
    # Pushforward (image measure)
    # ------------------------------------------------------------------

    def pushforward(
        self,
        phi: Callable[[np.ndarray], np.ndarray],
        *,
        new_space: Manifold,
    ) -> DiscreteMeasure:
        r"""Image measure :math:`\varphi_* \mu` under :math:`\varphi`.

        For any test function :math:`g` on the target space,

        .. math::
           :label: discrete-measure-pushforward

           \int_{\mathcal{Y}} g \, d(\varphi_* \mu)
           \;=\; \int_{\mathcal{X}} (g \circ \varphi) \, d\mu
           \;=\; \sum_i w_i \, g(\varphi(x_i)).

        Computationally this is "apply ``φ`` to the nodes, keep the
        weights." Note the asymmetry: the **change-of-variables
        Jacobian is NOT applied** — the pushforward is the φ-image
        of the discrete measure, not a Radon-Nikodym derivative
        against a target reference measure. If the user wants
        ``g(y) dy = h(x) dx`` semantics on a smooth target measure,
        they must pre-multiply the weights by the Jacobian
        :math:`|\det D\varphi|^{-1}` themselves.

        For a non-invertible ``φ``, :math:`\varphi_* \mu` may have
        atoms collapsing onto the same target point — the caller's
        ``g`` will then see those weights summed implicitly through
        the integration step. Mathematically valid; numerically the
        node array will contain duplicates with separate weights.

        Parameters
        ----------
        phi : callable
            Map :math:`\varphi : \mathcal{X} \to \mathcal{Y}`.
            Called once with the full ``nodes`` array. Must return
            an array of shape ``(N,)`` or ``(N, d')``.
        new_space : Manifold
            The target point set :math:`\mathcal{Y}` — **required**.

            ⭐ It was optional until 2026-09-01, defaulting to a
            *fabricated* support named ``f"φ_*({self.support})"``. That
            default is the defect this campaign exists to remove: the image
            of a manifold under an arbitrary map is a manifold nobody has
            computed, and naming it does not make it known. Only ``φ``'s
            author knows :math:`\mathcal{Y}`, so only they can name it.
            (``[M]`` at the retype: 7 of 8 call sites already passed it; the
            one that did not was a test asserting the fabricated name.)

        Returns
        -------
        DiscreteMeasure
            A new measure with ``nodes = φ(nodes)`` and unchanged
            weights. ``invariance_group`` is dropped (in general
            ``φ`` does not preserve a group action), and
            ``degree_of_exactness`` is dropped unless ``φ`` is the
            identity (which we cannot statically determine —
            caller's responsibility to set if desired).
        """
        new_nodes = np.asarray(phi(self.nodes))
        if new_nodes.shape[0] != self.n_points:
            raise ValueError(
                f"pushforward map must preserve number of atoms, "
                f"got {new_nodes.shape[0]} from N={self.n_points}"
            )

        return DiscreteMeasure(
            nodes=new_nodes,
            weights=self.weights.copy(),
            support=new_space,
            invariance_group=None,
            exactness=None,
        )

    def consolidate(self, *, atol: float = 1e-12) -> DiscreteMeasure:
        r"""Merge coincident atoms, summing their weights.

        The reduction :meth:`pushforward` documents but does not perform:
        a non-invertible :math:`\varphi` collapses atoms onto a shared
        target point, leaving duplicate nodes carrying separate weights.
        That is mathematically valid — a caller integrating against it sums
        them implicitly — but it is an *unreduced* representation, and code
        that inspects ``nodes`` (counting them, partitioning them, matching
        them under a group action) sees a different object than the measure
        it represents.

        This is the second half of a quotient. With
        :math:`\varphi = ` "map each node to its orbit representative",

        .. math::

           \mathrm{quotient}(G) \;=\;
           \mathrm{pushforward}(\varphi)\ \texttt{.consolidate()},

        and the summed weights are exactly the orbifold measure
        :math:`W = w \cdot |G| / |\mathrm{Stab}|` that orbit-stabilizer
        predicts — derived, not chosen.

        **Total mass is preserved exactly** (the weights are summed, not
        dropped), which is the in-tree discriminator between a QUOTIENT and
        a :meth:`restrict`-style restriction: the latter drops mass.

        Unlike :meth:`pushforward`, both ``invariance_group`` and
        ``degree_of_exactness`` are **preserved**. Consolidation moves no
        node and changes no integral — :math:`\int g\,d\mu` is identical
        for every test function — so neither claim can be invalidated by
        it. A group that permuted the original atoms still permutes the
        merged ones, since merging is by position.

        Parameters
        ----------
        atol : float, optional
            Nodes within this distance are treated as one atom.

        Returns
        -------
        DiscreteMeasure
            The reduced measure. Idempotent: consolidating twice equals
            consolidating once.
        """
        nodes = self.nodes
        flat = nodes.reshape(self.n_points, -1)
        keys = (np.round(flat / atol) * atol + 0.0).astype(np.float64)

        # Group by rounded position, preserving first-appearance order so
        # the result is deterministic and reads like the input.
        seen: dict[bytes, int] = {}
        order: list[int] = []
        merged_weights: list[float] = []
        for i in range(self.n_points):
            key = keys[i].tobytes()
            slot = seen.get(key)
            if slot is None:
                seen[key] = len(order)
                order.append(i)
                merged_weights.append(float(self.weights[i]))
            else:
                merged_weights[slot] += float(self.weights[i])

        if len(order) == self.n_points:
            return self  # already reduced — nothing coincides

        return DiscreteMeasure(
            nodes=nodes[np.array(order, dtype=np.int64)],
            weights=np.array(merged_weights, dtype=np.float64),
            support=self.support,
            invariance_group=self.invariance_group,
            exactness=self.exactness,
            # Consolidation merges coincident atoms and changes no
            # integral, so the reference measure is untouched — the
            # same reason the degree survives.
        )

    # ------------------------------------------------------------------
    # Quotient by a symmetry group
    # ------------------------------------------------------------------

    def on_orbit_space(self, orbit_space: Quotient) -> DiscreteMeasure:
        r"""The same atoms, read as chart coordinates of an orbit space.

        A rule built on a manifold :math:`C` is a rule on every orbit
        space :math:`M/H` whose chart codomain
        (:attr:`~orpheus.numerics.manifold.Quotient.realization`) is
        :math:`C`: same points, same weights — only what the measure
        KNOWS about its support changes, from "an interval" to "the polar
        marginal of a sphere, about this axis". It is how the slab says
        its :math:`\mu`-rule lives on :math:`S^2/SO(2)_x` rather than on
        the interval a chart happens to map onto, which is ERR-080's
        defect class stated at the level of spaces (tracker 2.4, 2026-09-01):
        `[M]` before it, an 8-node slab ANGULAR space and an 8-node SPATIAL
        space on :math:`[-1,1]` were ``==`` and hash-equal.

        It is NOT a :meth:`pushforward` — no map is applied — and NOT a
        :meth:`quotient` — the fold starts from a measure on the BASE and
        identifies orbits, whereas a :math:`\mu`-rule was never on
        :math:`S^2` to begin with.

        ``invariance_group`` is DROPPED: a subgroup of :math:`O(3)` is a
        statement about an embedding, and the orbit space fixes an
        embedding (its axis) that the chart did not. The adopter that
        knows the residual re-tags it. ``exactness`` survives: the
        reference measure is a measure on the chart, and the chart is
        unchanged.

        Raises
        ------
        ValueError
            If ``orbit_space.realization`` is not this measure's support —
            the one precondition, refused where the declaration is written.
        """
        if orbit_space.realization != self.support:
            raise ValueError(
                f"a measure on {self.support.name!r} cannot be read on "
                f"{orbit_space.name!r}: that orbit space's chart is "
                f"{orbit_space.realization.name!r}. A rule declares the "
                f"orbit space whose CHART it was built on; to fold a rule on "
                f"the base manifold, use quotient()."
            )
        return replace(self, support=orbit_space, invariance_group=None)

    def quotient(
        self,
        group: "SubgroupOfO3",
        *,
        atol: float = 1e-13,
    ) -> DiscreteMeasure:
        r"""The quotient measure :math:`\mu / G` on the orbifold
        :math:`\mathcal{X}/G` — THE FOLD.

        Names the composite both halves already document (the equation
        renders as :eq:`discrete-measure-quotient` on the theory page —
        this module is not rendered by Sphinx, so the label lives
        there):

        .. math::

           \mu / G \;=\; \mathrm{pushforward}(\varphi_{\mathrm{rep}})
           \texttt{.consolidate()},
           \qquad
           \varphi_{\mathrm{rep}}(x_i) = x_{\min(G \cdot i)},

        where :math:`\varphi_{\mathrm{rep}}` sends every node to its
        orbit's representative. The representative is the orbit's
        first-appearing member — the only *group-generic* section (a
        geometric section such as :math:`\xi \mapsto |\xi|` exists only
        for a mirror; a rotation orbit has no sign to take the absolute
        value of) — so the quotient's atoms appear in the parent's own
        storage order.

        **The defining law** (what makes this a quotient and not a
        convention): for every :math:`G`-invariant :math:`f`,

        .. math::

           \int f \, d(\mu/G) \;=\; \int f \, d\mu ,

        because each orbit's weights are *summed onto* the
        representative, giving the orbit-stabilizer weight

        .. math::

           W \;=\; w \cdot |G| \, / \, |\mathrm{Stab}(x)|

        — derived, not chosen. For a non-invariant :math:`f` the
        quotient instead integrates :math:`f \circ \varphi_{\mathrm{rep}}`:
        a :math:`\xi`-odd integrand that vanished on the full sphere by
        cancellation is *reported differently* on the fold. That is the
        quotient being honest about its smaller space, not an error —
        consumers that need odd moments must keep the full measure.

        **Total mass is preserved exactly** (weights are summed, never
        dropped) — the in-tree discriminator between a QUOTIENT and a
        :meth:`restrict`-style restriction. When the action is **free**
        (empty singular set :math:`\Sigma = \varnothing`, i.e. no node
        is fixed by a non-identity element), every orbit has length
        exactly :math:`|G|` and the weights collapse to a uniform
        :math:`|G| \cdot w` — no mixed :math:`w`/:math:`2w` sum. A free
        action is the fold's well-posedness condition: a rule is
        admissible for a fold iff it places no node on :math:`\Sigma`.

        **Precondition by construction.** The quotient is defined only
        on a :math:`G`-invariant measure; the certificate
        (:func:`~orpheus.numerics.symmetry.orbit_certificate`) is the
        proof, and it cannot be built for a non-invariant measure, a
        CONTINUOUS group (no finite node permutation), or nodes without
        a 3-D realization — all three refuse loudly here.

        **The fold consumes the symmetry.** The result carries
        ``invariance_group=None``: the representative section is not
        :math:`G`-equivariant, so no residual invariance is guaranteed
        (ask :func:`~orpheus.numerics.symmetry.maximal_invariance_groups`
        of the result if one is needed). Consequently a SECOND
        ``quotient(group)`` refuses — except in the degenerate case
        where every node was fixed (the action was trivial), where the
        quotient is the identity on ``(nodes, weights)`` and literally
        idempotent. ``exactness`` is likewise dropped (through
        :meth:`pushforward`): an honest claim would be against the
        pushforward reference :math:`\varphi_* \lambda`, a type nothing
        consumes yet — minting it on speculation is worse than silence
        (the direct-sum precedent).

        Parameters
        ----------
        group : SubgroupOfO3
            The finite symmetry group to quotient by, e.g.
            ``SubgroupOfO3.Mirror("y")`` for the cylindrical
            :math:`\xi`-fold.
        atol : float, optional
            Node-matching tolerance for building the certificate — the
            one honestly-numerical step. Orbit membership itself is an
            integer identity on the resulting permutations.

        Returns
        -------
        DiscreteMeasure
            One atom per :math:`G`-orbit, at the orbit's first-appearing
            node, carrying the summed orbit weight; ``support`` tagged
            ``f"{support}/{group.name}"``.

        Raises
        ------
        ValueError
            If no certificate exists (measure not ``group``-invariant,
            or ``group`` continuous / without a finite realization on
            these nodes).
        """
        # Lazy import: ``symmetry`` imports FROM this module (it takes a
        # DiscreteMeasure), so a module-level import here would be
        # circular. Mirrors the ``space`` property's lazy import.
        from orpheus.numerics.symmetry import orbit_certificate

        certificate = orbit_certificate(self, group, atol=atol)
        if certificate is None:
            raise ValueError(
                f"a quotient is defined only for a {group.name}-invariant "
                f"measure with a finite realization; this measure is not "
                f"{group.name}-invariant (or {group.name} is continuous, "
                f"and a continuous group has no finite node permutation)"
            )

        representative = np.arange(self.n_points, dtype=np.int64)
        for orbit in certificate.orbits():
            representative[orbit] = orbit.min()

        return self.pushforward(
            lambda nodes: nodes[representative],
            # ⭐ The folded measure lives on the ORBIT SPACE, and the orbit space
            # is what ``Manifold.quotient`` builds — the same catalogue entry a
            # geometry consults, not a second spelling of its name.
            new_space=self.support.quotient(group),
        ).consolidate()

    # ------------------------------------------------------------------
    # Restriction
    # ------------------------------------------------------------------

    def restrict(
        self,
        predicate: Callable[[np.ndarray], np.ndarray],
    ) -> DiscreteMeasure:
        r"""Restriction :math:`\mathbf{1}_E \cdot \mu` to the set
        :math:`E = \{ x : \text{predicate}(x) \text{ is true} \}`.

        Drops atoms outside :math:`E`; weights of kept atoms are
        preserved. Used for half-range SN sweeps
        (:math:`E = \{\mu > 0\}`), for cutting MoC bundles to a
        region, and generally for any measure restricted to a
        measurable subset.

        Parameters
        ----------
        predicate : callable
            Vectorised mask function. Called once with the full
            ``nodes`` array; must return a boolean array of shape
            ``(N,)``.

        Returns
        -------
        DiscreteMeasure
            A new measure with the kept subset of nodes/weights.
            ``space`` is unchanged (the restriction lives on the
            same space, with support a subset). ``invariance_group``
            is dropped (a restriction typically breaks any invariance
            the original carried, unless ``E`` is invariant — caller
            may re-tag).  ``degree_of_exactness`` is dropped (the
            restricted rule is no longer exact for the original
            polynomial class).
        """
        mask = np.asarray(predicate(self.nodes), dtype=bool)
        if mask.shape != (self.n_points,):
            raise ValueError(
                f"predicate must return a boolean array of shape "
                f"({self.n_points},), got shape {mask.shape}"
            )
        return DiscreteMeasure(
            nodes=self.nodes[mask],
            weights=self.weights[mask],
            support=self.support,
            invariance_group=None,
            exactness=None,
        )

    # ------------------------------------------------------------------
    # Partition (disjoint decomposition by labelling predicate)
    # ------------------------------------------------------------------

    def partition_by(
        self,
        predicate: Callable[[np.ndarray], np.ndarray],
    ) -> tuple["DiscreteMeasurePartition", ...]:
        r"""Disjoint partition by node-label predicate.

        Realises the inverse of :meth:`__add__`: from a measure
        :math:`\mu` and a labelling map :math:`\ell : \mathcal{X} \to L`
        on the nodes, returns the disjoint decomposition

        .. math::
           :label: discrete-measure-partition

           \mu \;=\; \bigoplus_{\lambda \in L} \mu_\lambda,
           \qquad
           \mu_\lambda
           \;=\; \sum_{i : \ell(x_i) = \lambda} w_i \, \delta_{x_i}.

        Each output entry carries the partition label, the indices
        into the parent measure, and the restricted measure
        :math:`\mu_\lambda`. The decomposition is **disjoint** by
        construction (every node appears in exactly one partition
        entry) and **preserves total mass** (sum of partition weights
        equals total mass of :math:`\mu`).

        The canonical use is the **angular octant partition** of a
        spherical cubature: the predicate
        :math:`\ell(\Omega) = (\mathrm{sign}\,\mu_x, \mathrm{sign}\,\mu_y,
        \mathrm{sign}\,\mu_z)` produces the eight octants of :math:`S^2`
        (or four for a 2-D cubature where :math:`\mu_z = 0` is degenerate
        and contributes a separate ``sign=0`` entry).

        Other consumers (per :ref:`tensorial-framing`):

        * MoC track-bundle direction grouping (per polar angle).
        * MC boundary-current scoring by hemisphere.
        * SN boundary-realiser incoming-ordinate masking
          (Grand Report v3 §16A.5).

        Parameters
        ----------
        predicate : callable
            Vectorised labelling map. Called once with the full
            ``nodes`` array (shape ``(N,)`` for 1-D nodes or
            ``(N, d)`` for multi-dimensional). Must return either:

            - A 1-D array of shape ``(N,)`` of scalar labels (e.g.
              integer octant indices), OR
            - A 2-D array of shape ``(N, k)`` of compound labels
              (e.g. octant signs ``(sign_x, sign_y, sign_z)``).

        Returns
        -------
        tuple of DiscreteMeasurePartition
            One entry per unique label, in lexicographic order of
            labels. Each entry carries:

            - ``label`` — the predicate output, as a tuple (always a
              tuple even for scalar labels, so keys round-trip cleanly).
            - ``indices`` — 1-D ``int`` array of indices into
              ``self.nodes`` / ``self.weights``.
            - ``measure`` — the restricted :class:`DiscreteMeasure`
              on the same space as ``self``.

        Notes
        -----

        Inverse-of-direct-sum identity: for any partition predicate
        :math:`\ell`, the direct sum of partition measures equals
        ``self`` modulo ordering and metadata,

        .. math::

           \mu \;=\; \mu_{\lambda_1} \oplus \mu_{\lambda_2}
                     \oplus \cdots \oplus \mu_{\lambda_K}.

        ``invariance_group`` and ``degree_of_exactness`` are dropped
        on each partition entry — a partition typically breaks any
        global invariance, and the polynomial-exactness of the parent
        rule does not survive restriction.
        """
        labels = np.asarray(predicate(self.nodes))
        if labels.shape[0] != self.n_points:
            raise ValueError(
                f"predicate must return labels for all {self.n_points} "
                f"nodes; got leading dimension {labels.shape[0]}."
            )
        if labels.ndim == 1:
            unique_labels = np.unique(labels)
            label_tuples = [(_to_python_scalar(u),) for u in unique_labels]
            masks = [labels == u for u in unique_labels]
        elif labels.ndim == 2:
            # Compound label: group by row equality.
            unique_rows = np.unique(labels, axis=0)
            label_tuples = [
                tuple(_to_python_scalar(v) for v in row) for row in unique_rows
            ]
            masks = [np.all(labels == row, axis=1) for row in unique_rows]
        else:
            raise ValueError(
                f"predicate output must be 1-D or 2-D, got shape "
                f"{labels.shape}"
            )

        result = []
        for label, mask in zip(label_tuples, masks):
            indices = np.where(mask)[0]
            measure = DiscreteMeasure(
                nodes=self.nodes[indices],
                weights=self.weights[indices],
                support=self.support,
                invariance_group=None,
                exactness=None,
            )
            result.append(
                DiscreteMeasurePartition(
                    label=label, indices=indices, measure=measure,
                )
            )
        return tuple(result)

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def with_metadata(
        self,
        *,
        invariance_group: "SubgroupOfO3 | None" = None,
        exactness: "ExactnessClaim | None" = None,
    ) -> DiscreteMeasure:
        """Return a copy with optional metadata fields populated.

        Used by Issue 3 (symmetry tagging) and by quadrature-rule
        constructors that tag a measure after building it.

        The former ``degree_of_exactness: int`` parameter is **retired**,
        not renamed: it took a bare degree with no reference, so it was
        an API through which a half-claim could be attached to a rule
        after the fact — the exact state the exactness carve exists to
        make unspellable. `[M]` it had **zero** callers; every
        ``with_metadata`` call in the tree passes ``invariance_group``
        only.
        """
        return replace(
            self,
            invariance_group=(
                invariance_group
                if invariance_group is not None
                else self.invariance_group
            ),
            exactness=(
                exactness if exactness is not None else self.exactness
            ),
        )


# ---------------------------------------------------------------------------
# Partition entry (return type of DiscreteMeasure.partition_by)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DiscreteMeasurePartition:
    r"""One entry of a :meth:`DiscreteMeasure.partition_by` decomposition.

    Carries a partition label, the indices into the parent measure,
    and the restricted measure :math:`\mu_\lambda` on the same space.

    The label is a Python tuple — even for scalar labels it is wrapped
    as a length-1 tuple so partition-keyed dicts (e.g.
    ``Mapping[OctantLabel, SweepDependencyGraph]`` in the SN sweep)
    have uniform key shape regardless of label dimensionality.

    Attributes
    ----------
    label : tuple
        The predicate output for this partition. For a scalar-label
        predicate (e.g. integer index), a length-1 tuple
        ``(idx,)``. For a compound-label predicate (e.g. octant signs
        ``(sign_x, sign_y, sign_z)``), the full tuple.
    indices : np.ndarray
        Indices into the parent measure's ``nodes`` / ``weights``.
        Shape ``(N_part,)``, dtype ``int``.
    measure : DiscreteMeasure
        Restricted measure :math:`\mu_\lambda` on the same space as
        the parent. ``invariance_group`` and ``degree_of_exactness``
        are dropped (a partition typically breaks any invariance the
        parent had).
    """

    label: tuple
    indices: np.ndarray
    measure: DiscreteMeasure


def _to_python_scalar(v):
    """Convert a numpy scalar to a built-in int/float for tuple labels.

    Used internally by :meth:`DiscreteMeasure.partition_by` so the
    returned ``label`` tuples have stable Python types (suitable for
    dict keys).
    """
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return float(v)
    return v


# ---------------------------------------------------------------------------
# Bundle measures (disintegration)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BundleMeasure:
    r"""Disintegrated measure on a fibered space.

    A fibered space is a measurable space :math:`\mathcal{X}` equipped
    with a projection :math:`\pi : \mathcal{X} \to \mathcal{B}` to a
    base space :math:`\mathcal{B}`. By the disintegration theorem
    (Bourbaki 1969, Integration VI §3), under mild regularity any
    measure :math:`\mu` on :math:`\mathcal{X}` decomposes as

    .. math::
       :label: bundle-measure-disintegration

       \int_\mathcal{X} f \, d\mu
       \;=\; \int_\mathcal{B} \! \left[
           \int_{\pi^{-1}(b)} f \, d\mu_b
       \right] d\mu_\mathcal{B}(b),

    i.e. a base measure :math:`\mu_\mathcal{B}` on :math:`\mathcal{B}`
    plus, for each base point :math:`b`, a fiber measure :math:`\mu_b`
    on the fiber :math:`\pi^{-1}(b)`.

    For MoC ray-bundle quadratures (Wave 2 of the SN-reshape campaign)
    the base is the angular sphere (or a half-range thereof) and the
    fiber over each direction is the set of parallel rays through the
    geometry — a set whose discretisation depends on the direction.
    The fiber measures cannot be hoisted into a single product
    measure because their support varies with the base point.

    SN does not consume :class:`BundleMeasure` directly in the SN
    reshape campaign, but the abstraction is shipped in Phase 0 so
    MoC migration in Wave 2 does not have to revisit the primitives
    layer (per the design note in
    ``.claude/plans/sn_reshape.md`` Issue 2).

    Parameters
    ----------
    base : DiscreteMeasure
        Measure on the base space :math:`\mathcal{B}`.
    fiber_factory : callable
        Map sending each base node ``b`` (a single row of
        ``base.nodes``) to the fiber measure :math:`\mu_b` on
        :math:`\pi^{-1}(b)`. Called lazily — only when
        :meth:`integrate` is invoked or when explicit fiber
        materialisation is requested. Returns a
        :class:`DiscreteMeasure`.

    Notes
    -----

    The class is **frozen**. Composition operations (tensor product
    of bundles, concatenation across bases, …) are not implemented
    here — they are not needed by SN or MoC at Phase 0. Add them
    incrementally when a consumer needs them.
    """

    base: DiscreteMeasure
    fiber_factory: Callable[[np.ndarray], DiscreteMeasure]

    def integrate(
        self,
        f: Callable[[np.ndarray, np.ndarray], np.ndarray],
    ) -> np.ndarray:
        r"""Iterated integral against the disintegrated measure.

        Computes :math:`\int_\mathcal{B} \int_{\pi^{-1}(b)} f(b, x) \,
        d\mu_b(x) \, d\mu_\mathcal{B}(b)`. See
        :eq:`bundle-measure-disintegration`.

        Parameters
        ----------
        f : callable
            Two-argument test function ``f(base_node, fiber_node)``
            returning a scalar (or array). The base node is a single
            row of ``base.nodes`` (scalar if ``base.dim == 1``); the
            fiber node argument is the **full fiber-node array** for
            that base point — the function is expected to be
            vectorised over the fiber axis.

        Returns
        -------
        np.ndarray
            The iterated sum :math:`\sum_i w_i^B \sum_j w_{ij}^F
            f(b_i, x_{ij})`.
        """
        base_nodes = self.base.nodes
        is_1d = base_nodes.ndim == 1

        total: float | np.ndarray = 0.0
        for i, w_b in enumerate(self.base.weights):
            b = base_nodes[i] if is_1d else base_nodes[i]
            fiber = self.fiber_factory(b)
            inner = fiber.integrate(lambda x, b=b: f(b, x))
            total = total + w_b * inner
        return np.asarray(total)


# ---------------------------------------------------------------------------
# 1-D primitive constructors
# ---------------------------------------------------------------------------


def gauss_legendre(n: int) -> DiscreteMeasure:
    r"""N-point Gauss-Legendre quadrature on :math:`[-1, 1]`.

    Returns the discrete measure :math:`\mu = \sum_{i=1}^n w_i \,
    \delta_{x_i}` whose nodes :math:`x_i` are the roots of the
    Legendre polynomial :math:`P_n` and whose weights are the
    Christoffel numbers
    :math:`w_i = 2 \big/ \bigl[(1 - x_i^2) (P_n'(x_i))^2\bigr]`.

    Polynomial exactness (Stoer & Bulirsch 2002, Theorem 3.6.20):
    exact for every polynomial of degree :math:`\le 2n - 1`. The
    returned measure carries ``degree_of_exactness = 2*n - 1`` so
    downstream consumers can verify the claim.

    Weight sum: :math:`\sum_i w_i = 2`.

    Parameters
    ----------
    n : int
        Number of quadrature nodes; must be :math:`\ge 1`.

    Returns
    -------
    DiscreteMeasure
        On ``support=COSINE_INTERVAL`` with ``degree_of_exactness=2n-1``.

    See Also
    --------
    :meth:`orpheus.numerics.quadrature.Quadrature.gauss_legendre` — the
    named factory SN consumers call. There is no per-family adapter
    class: the four SN-side wrappers this docstring used to point at
    (``orpheus.sn.quadrature.GaussLegendre1D`` and its siblings) were
    retired into classmethod factories on the one ``Quadrature`` type.
    """
    # Late import: ``generating_measure`` imports DiscreteMeasure from
    # this module, so the runtime arrow has to run one way. See the
    # TYPE_CHECKING block at the top.
    from orpheus.numerics.generating_measure import LEGENDRE

    return LEGENDRE.gauss(n)


def gauss_chebyshev(n: int) -> DiscreteMeasure:
    r"""N-point Gauss-Chebyshev (first-kind) quadrature on
    :math:`[-1, 1]` with weight :math:`(1 - x^2)^{-1/2}`.

    The closed-form rule (Stoer & Bulirsch 2002, §3.6) is

    .. math::
       :label: discrete-measure-gauss-chebyshev

       x_i = \cos\!\left( \frac{(2i - 1) \pi}{2n} \right),
       \qquad w_i = \frac{\pi}{n}, \qquad i = 1, \ldots, n,

    and the rule is exact in the weighted sense

    .. math::

       \sum_{i=1}^n w_i \, q(x_i)
       \;=\; \int_{-1}^{1} \frac{q(x)}{\sqrt{1 - x^2}} \, dx

    for every polynomial :math:`q` of degree :math:`\le 2n - 1`.

    Note the **weight function is in the integral**, not in the
    quadrature — :math:`\int q \, d\mu_{\text{GC}} = \int q(x) (1 -
    x^2)^{-1/2} dx`, so an unweighted polynomial integration on
    :math:`[-1, 1]` does **not** simplify to evaluating
    :math:`\int q \, dx`. Callers wanting unweighted polynomial
    integration should use :func:`gauss_legendre`.

    Weight sum: :math:`\sum_i w_i = \pi`.

    Parameters
    ----------
    n : int
        Number of quadrature nodes; must be :math:`\ge 1`.

    Returns
    -------
    DiscreteMeasure
        On ``support=COSINE_INTERVAL`` with ``degree_of_exactness=2n-1``
        (with respect to the weighted integral).
    """
    from orpheus.numerics.generating_measure import CHEBYSHEV_T

    return CHEBYSHEV_T.gauss(n)


def equispaced(a: float, b: float, n: int) -> DiscreteMeasure:
    r"""Equispaced (midpoint-rule) measure on :math:`[a, b]`.

    Returns the discrete measure with :math:`n` nodes at the midpoints
    of an equal partition of :math:`[a, b]`, all carrying weight
    :math:`(b - a) / n`:

    .. math::

       x_i = a + \left(i - \tfrac{1}{2}\right)\,\frac{b - a}{n},
       \qquad w_i = \frac{b - a}{n}.

    The midpoint rule has degree of exactness ``1`` (exact for
    constants and linears, not for quadratics in general).

    This primitive is offered for tensor-product azimuthal quadrature,
    but it is **not** what the shipped SN product rule uses. The φ-axis
    of :func:`~orpheus.numerics.quadrature.product_mu_phi` is
    :func:`~orpheus.numerics.quadrature.periodic_trapezoid` at
    ``shift=NODE_ALIGNED`` — the left-endpoint grid
    :math:`\varphi_m = 2\pi m / n_\varphi`, generated as roots of
    unity, not midpoints. (That convention was carried by the SN-side
    ``ProductQuadrature`` adapter this docstring used to point at;
    R-1 Phase A detour-C retired all four per-family adapters into
    classmethod factories on the one ``Quadrature`` type, and the
    convention now lives in the rule function.) Midpoints are offered
    here because they integrate constants exactly while preserving
    symmetry under reflection through the centre of the interval —
    but the shift a rule *needs* is a well-posedness question (whether
    a node sits on the mirror axis), so select it deliberately rather
    than reaching for this primitive by default.

    Parameters
    ----------
    a, b : float
        Endpoints of the interval; ``a < b`` is required.
    n : int
        Number of nodes; must be :math:`\ge 1`.

    Returns
    -------
    DiscreteMeasure
        On ``support=Interval(a, b)`` with ``degree_of_exactness=1``.
    """
    if n < 1:
        raise ValueError(f"equispaced requires n >= 1, got n={n}")
    if not (a < b):
        raise ValueError(f"equispaced requires a < b, got a={a}, b={b}")
    from orpheus.numerics.exactness import (
        ExactnessClaim,
        OrthogonalSystem,
        UniformMeasure,
    )

    h = (b - a) / n
    nodes = a + (np.arange(n) + 0.5) * h
    weights = np.full(n, h)
    support = Interval(a, b)
    return DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=support,
        # ALGEBRAIC degree 1, against Lebesgue measure on THIS INTERVAL —
        # and naming both halves is the point. The same nodes read as a
        # rule on the CIRCLE are the periodic trapezoid, exact to
        # trigonometric degree n-1; that is a different claim about a
        # different reference, and it belongs to a circle rule rather
        # than to this one (§4: a rule on a domain is ONE object).
        exactness=ExactnessClaim(
            reference=UniformMeasure(
                support=support,
                orthogonal_system=OrthogonalSystem.ALGEBRAIC,
            ),
            degree=1,
        ),
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


__all__ = [
    "BundleMeasure",
    "DiscreteMeasure",
    "DiscreteMeasurePartition",
    "equispaced",
    "gauss_chebyshev",
    "gauss_legendre",
]

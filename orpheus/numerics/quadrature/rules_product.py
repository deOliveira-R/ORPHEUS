r"""Product quadrature on the sphere: Gauss-Legendre :math:`(\mu)`
:math:`\times` equispaced :math:`(\phi)`.

Both factors are **registered rules**, not inline constructions: the
polar one is :func:`~orpheus.numerics.quadrature.rules_1d.gauss_legendre_on_mu`
on :math:`\mu = \cos\theta \in [-1, 1]`, and the azimuthal one is
:func:`~orpheus.numerics.quadrature.rules_circle.periodic_trapezoid` on
:math:`S^1` at ``shift=NODE_ALIGNED`` — the :math:`n_\phi`
left-endpoints of an even partition of :math:`[0, 2\pi)`, which is the
long-standing convention here, now *named* rather than spelled as a
``linspace`` literal.

That naming is the seam a fold needs: selecting the **staggered** shift
(and with it :math:`\Sigma = \emptyset`) becomes a substitution of one
registered rule for another, rather than a new flag on this function.
It also makes both factors carry their own exactness claims, which is
what lets this rule's degree be *derived* — see
:func:`spherical_product_claim`.

Direction cosines per ordinate:

.. math::
   :label: product-mu-phi-cosines

   \mu_x &= \sin\theta \cos\phi, \\
   \mu_y &= \sin\theta \sin\phi, \\
   \mu_z &= \cos\theta = \mu.

Weight per ordinate: :math:`w = w_{\text{GL}}(\mu) \cdot (2\pi /
n_\phi)`. Total weight sum: :math:`4\pi`.

Symmetry: :math:`D_{n_\phi h}`, of order :math:`4 n_\phi`. The
:math:`n_\phi` equispaced azimuths are preserved by the cyclic
rotations :math:`C_{n_\phi}` and by the vertical mirrors at
:math:`k\pi/n_\phi`; the Gauss-Legendre :math:`\mu` nodes are
symmetric about :math:`\mu = 0`, which supplies :math:`\sigma_h`.
Together these generate :math:`D_{n_\phi h}` — and the tag is now
**derived at construction** from exactly those three generator facts,
each measured on its factor through
:meth:`~orpheus.geometry.transformation.RigidMotion.preserves` (see
:func:`spherical_product`); a factor pair failing a check is refused,
never mis-tagged. The independent realization —
:func:`~orpheus.numerics.symmetry.maximal_invariance_groups`
computing the group from the ASSEMBLED nodes — pins the same answer
in ``tests/``, so the derivation and the walk verify each other.

.. caution::

   The agreement holds to ~1e-16, not bit-exactly, and the match window
   (``symmetry._NODE_WINDOW_FACTOR``) is what absorbs the difference.
   Lebedev's :math:`O_h` claim, by contrast, is exact on both sides
   because signed permutations are exact in IEEE.

   Issue #325 named **two** halves for closing this gap, and BOTH are
   done: the :math:`\phi` grid is
   :mod:`~orpheus.numerics.roots_of_unity`-generated (2026-08-02), and
   the **checker's** own :math:`C_n` operators consume the SAME
   generator (``symmetry._cyclic_ops`` builds each rotation from an
   exact circle point, and the :math:`\sigma_v` mirrors are the coset
   :math:`C_n \cdot \sigma_0` of a bit-exact signed diagonal — one
   trigonometric spelling in the whole chain, not two). The window
   survives for what remains genuinely numerical: the matrix-times-node
   MULTIPLICATION round-off when an operator is applied to the nodes.

.. note::

   Until 2026-08-02 this rule tagged :math:`SO(2)` and this
   paragraph defended it as "a conservative upper bound". Both
   halves were wrong, in opposite directions. **No finite point set
   on** :math:`S^2` **is** :math:`SO(2)`-**closed**, so the tag was
   simply false — and for an *invariance* claim a larger group is a
   **stronger** claim, so an upper bound is an over-claim, never a
   conservative one. (Contrast ``degree_of_exactness``, where the
   ``min`` genuinely *is* conservative: a lower degree is a weaker
   claim.) The prose also undersold the truth in the other
   direction, calling the rule ":math:`C_{n_\phi}`-invariant
   strictly" — :math:`C_{n_\phi}` is an index-4 subgroup of the real
   answer. The false tag was load-bearing: it was the only reason
   the registry's geometry gate admitted this rule for a cylinder
   (see :func:`~orpheus.numerics.quadrature.registry.select_quadrature`).

References
----------

* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010). "The Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes."
  *Nucl. Sci. Eng.* **165**(2), 149-169, doi:10.13182/NSE08-66.
  **Eq. (50)** is the cylindrical (R-Z) :math:`\alpha`-recursion used by
  the sweep that consumes this product rule; **Eq. (52)** fixes the
  order of a level, which is what this module's per-level index lists
  are sorted by.

  ⛔ *Retracted citation (2026-08-27).* This entry named a
  non-existent "Bailey, Adams, Yang, Zika (2009), Annals of Nuclear
  Energy 35, 1929-1936"; the equation numbers were right and belong to
  BMC 2010 above. Full account + the ⚠ published-typo warning on
  Eq. (50): ``docs/theory/methods/sn/curvilinear_one_group.rst`` §bmc-equation-map.
"""

from __future__ import annotations

import numpy as np

from orpheus.geometry.transformation import RigidMotion

from ..exactness import (
    UNIFORM_ON_SPHERE,
    ExactnessClaim,
    OrthogonalSystem,
)
from ..measure import SPACE_CIRCLE, SPACE_SPHERE, DiscreteMeasure
from ..symmetry import SubgroupOfO3
from .rules_1d import gauss_legendre_on_mu
from .rules_circle import NODE_ALIGNED, periodic_trapezoid
from .rules_sphere import LevelStructure, PolarInvariant


def spherical_product_claim(
    polar: ExactnessClaim, azimuthal: ExactnessClaim
) -> ExactnessClaim:
    r"""The claim of a polar :math:`\times` azimuthal product on :math:`S^2`.

    .. math::
       :label: spherical-product-degree

       \deg_{S^2} \;=\; \min\big(\deg_{\text{polar}},\,
                                 \deg_{\text{azimuthal}}\big),

    against the uniform measure on :math:`S^2`, in the
    **spherical-harmonic** system.

    Why this is not
    :meth:`~orpheus.numerics.exactness.ExactnessClaim.tensor_with`
    ---------------------------------------------------------------

    The factors' systems differ — algebraic in :math:`\mu`, trigonometric
    in :math:`\varphi` — so ``tensor_with`` correctly returns ``None``
    for this pair. It is right to: a tensor product lands on the
    **square** :math:`[-1,1] \times S^1`, spanning a mixed tensor system.
    This rule instead composes the factors through the **embedding**
    :math:`(\mu, \varphi) \mapsto \Omega`, so its claim is about
    spherical harmonics on :math:`S^2` — a different space, and therefore
    a different theorem. This is the theorem
    :class:`~orpheus.numerics.exactness.ProductMeasure`'s docstring means
    by "belongs with the rule that applies it".

    The derivation
    --------------

    A product rule **factorises** on a separated integrand, so it is
    enough to check the monomials :math:`\Omega_x^a \Omega_y^b
    \Omega_z^c` of total degree :math:`d`. Writing :math:`k = a + b` for
    the transverse power, such a monomial is
    :math:`(1-\mu^2)^{k/2}\,\mu^{\,d-k}` times a trigonometric polynomial
    of degree :math:`k` in :math:`\varphi`, and the rule's value is the
    product of the two factors' values. Two cases:

    * :math:`k` **odd** — the polar factor :math:`(1-\mu^2)^{k/2}` is
      **not a polynomial in** :math:`\mu`, so the polar rule has no
      exactness to offer. It does not need any: :math:`\cos^a\varphi
      \sin^b\varphi` with :math:`a + b` odd contains only *odd*
      harmonics, so every :math:`m \neq 0` and the azimuthal rule
      reproduces the exact :math:`0` — provided its degree is
      :math:`\ge k`. **The azimuthal factor annihilates exactly the terms
      whose polar factor is non-polynomial.**
    * :math:`k` **even** — the polar factor **is** a polynomial, of
      degree exactly :math:`d`, so the polar rule must be exact to
      :math:`d`; the azimuthal trigonometric degree is :math:`k \le d`.

    Both conditions read :math:`\deg \ge d`, so the largest :math:`d`
    both admit is the minimum — :eq:`spherical-product-degree`. `[M]` The
    bound is attained: ``product_mu_phi(4, 8)`` reproduces every
    degree-7 spherical monomial to **1.1e-17** and misses at degree 8 by
    **7.3e-2** (worst at :math:`\Omega_z^8`).

    That degree-7 residual was **1.1e-16** until 2026-08-02 — an order of
    magnitude larger — and the exact azimuths are what closed it. The
    tightness at degree 8 is unmoved, as it must be: it is a statement
    about which subspace the rule integrates, not about how well.

    Parameters
    ----------
    polar : ExactnessClaim
        The polar factor's own claim. MUST be
        :attr:`~orpheus.numerics.exactness.OrthogonalSystem.ALGEBRAIC` —
        the derivation's even-:math:`k` case is a statement about
        polynomials in :math:`\mu`.
    azimuthal : ExactnessClaim
        The azimuthal factor's own claim. MUST be
        :attr:`~orpheus.numerics.exactness.OrthogonalSystem.TRIGONOMETRIC`
        — the odd-:math:`k` case rests on the rule annihilating every
        non-zero harmonic, which is what a trigonometric degree asserts
        and an algebraic one does not.

    Returns
    -------
    ExactnessClaim
        Against :data:`~orpheus.numerics.exactness.UNIFORM_ON_SPHERE`, in
        the spherical-harmonic system.

    Raises
    ------
    ValueError
        If either factor is in the wrong orthogonal system. The guard is
        the point of typing the inputs: passing two *algebraic* claims
        (the shape a bare-integer ``min`` cannot distinguish) is the
        square's tensor product, not this embedding, and its degree is
        not this minimum.
    """
    if polar.system is not OrthogonalSystem.ALGEBRAIC:
        raise ValueError(
            f"the spherical product theorem needs an ALGEBRAIC polar "
            f"claim, got {polar.system.value} ({polar}). Its even-k case "
            f"is a statement about polynomials in mu."
        )
    if azimuthal.system is not OrthogonalSystem.TRIGONOMETRIC:
        raise ValueError(
            f"the spherical product theorem needs a TRIGONOMETRIC "
            f"azimuthal claim, got {azimuthal.system.value} "
            f"({azimuthal}). Its odd-k case rests on the rule "
            f"annihilating every non-zero harmonic — a trigonometric "
            f"degree asserts that; an algebraic one does not. The same "
            f"nodes read as an interval rule carry degree 1, which is "
            f"how a naive min() over bare integers reported 1 here."
        )
    return ExactnessClaim(
        reference=UNIFORM_ON_SPHERE,
        degree=min(polar.degree, azimuthal.degree),
    )


def _derived_product_group(
    polar: DiscreteMeasure,
    azimuthal: DiscreteMeasure,
    *,
    atol: float,
) -> SubgroupOfO3:
    r"""The product measure's invariance group, COMPUTED from the factors.

    Three generator checks through the verified
    :meth:`~orpheus.geometry.transformation.RigidMotion.preserves` core,
    each on the FACTOR it concerns (small arrays — cheap at
    construction, unlike an orbit closure over the assembled sphere
    rule):

    * the azimuthal node set is closed under rotation by
      :math:`2\pi/n_\varphi` — the :math:`C_{n_\varphi}` generator
      (closure under a generator implies closure under its powers, so a
      single check is sound here — this is NOT the ERR-072 trap, which
      sampled a CONTINUOUS group by a finite subset);
    * the azimuthal node set is closed under the mirror
      :math:`y \to -y` — the standard-setting vertical mirror, which
      together with :math:`C_{n_\varphi}` generates all
      :math:`n_\varphi` vertical mirrors at :math:`k\pi/n_\varphi`;
    * the polar node set is closed under :math:`\mu \to -\mu` — which
      the embedding turns into :math:`\sigma_h`.

    All three together generate :math:`D_{n_\varphi h}` (order
    :math:`4 n_\varphi`) in the standard setting. The checks are run on
    the factors rather than inherited from their ``invariance_group``
    tags: this module shipped three false symmetry declarations
    (ERR-072/073/074), so the derivation's inputs are *measured* facts,
    not declarations. The independent realization — the Hasse walk
    :func:`~orpheus.numerics.symmetry.maximal_invariance_groups` over
    the ASSEMBLED nodes — pins the same answer in ``tests/``; the two
    verify each other.

    Raises
    ------
    ValueError
        If the checks derive a group family ``SubgroupOfO3`` does not
        realize (:math:`C_{nv}` — no :math:`\sigma_h`;
        :math:`C_{nh}` — no vertical mirror; bare :math:`C_n`), or if
        the azimuthal factor is not even :math:`C_{n_\varphi}`-closed.
        Refusing beats tagging a wrong or weaker group: the field means
        MAXIMAL, so an under-claim is as false as an over-claim. No
        shipped factor pair reaches these arms.
    """
    n_phi = azimuthal.n_points
    e1 = np.array([1.0, 0.0])
    e2 = np.array([0.0, 1.0])
    step = RigidMotion.rotation(plane=(e1, e2), angle=2.0 * np.pi / n_phi)
    vertical_mirror = RigidMotion.reflection(normal=e2)
    polar_mirror = RigidMotion.reflection(normal=np.array([1.0]))

    rotation_ok = (
        step.preserves(
            azimuthal.nodes, azimuthal.weights, atol=atol, weight_atol=atol
        )
        is not None
    )
    if not rotation_ok:
        raise ValueError(
            f"the azimuthal factor is not invariant under rotation by "
            f"2π/{n_phi} — a spherical product's fiber must be a "
            f"C_{n_phi}-closed circle rule (equal weights on a closed "
            f"node set), got {azimuthal.n_points} nodes on "
            f"{azimuthal.support!r} that are not"
        )
    mirror_ok = (
        vertical_mirror.preserves(
            azimuthal.nodes, azimuthal.weights, atol=atol, weight_atol=atol
        )
        is not None
    )
    polar_ok = (
        polar_mirror.preserves(
            polar.nodes.reshape(-1, 1),
            polar.weights,
            atol=atol,
            weight_atol=atol,
        )
        is not None
    )
    if mirror_ok and polar_ok:
        return SubgroupOfO3.Dnh(n_phi)

    missing = (
        f"C_{n_phi}v" if mirror_ok else f"C_{n_phi}h" if polar_ok else f"C_{n_phi}"
    )
    raise ValueError(
        f"the product's invariance group would be {missing}, a family "
        f"SubgroupOfO3 does not realize (polar mirror-symmetric="
        f"{polar_ok}, azimuthal mirror-symmetric={mirror_ok}). Refusing "
        f"to tag rather than declare a wrong or weaker group — the "
        f"field means MAXIMAL. No shipped factor pair reaches this arm; "
        f"realize the family when a consumer does."
    )


def spherical_product(
    polar: DiscreteMeasure,
    azimuthal: DiscreteMeasure,
    *,
    atol: float = 1e-13,
) -> tuple[DiscreteMeasure, LevelStructure]:
    r"""Compose a polar rule and a circle rule into a rule on :math:`S^2`.

    The constructive twin of :func:`spherical_product_claim` — that
    function composes the factors' *claims* through the product theorem;
    this one composes their *measures* through the embedding
    :eq:`product-mu-phi-cosines`,

    .. math::

       (\mu, \varphi) \;\mapsto\; \Omega
       = (\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \mu),
       \qquad \sin\theta = \sqrt{1 - \mu^2},

    with product weights :math:`w = w_{\text{polar}} \cdot
    w_{\text{azimuthal}}`.

    **This signature is the seam the fold needs** (campaign Q5.2): the
    azimuthal factor arrives as an argument, so selecting the staggered
    circle rule — and with it an empty singular set
    :math:`\Sigma = \emptyset`, the fold's well-posedness condition —
    is a substitution of one registered rule for another::

        spherical_product(
            gauss_legendre_on_mu(n_mu),
            periodic_trapezoid(n_phi, shift=STAGGERED),
        )

    No ``half_range=True`` flag exists, because none is needed: the
    boolean the flag would encode is a property of the azimuthal
    *measure* (its shift classification), and it travels with the rule.

    Every intrinsic property of the result is DERIVED from the factors:

    * ``exactness`` through :func:`spherical_product_claim` (the
      theorem refuses mismatched systems — an interval rule handed in
      as the fiber is caught there, which is E3's same-nodes /
      wrong-space confusion made unrepresentable);
    * ``invariance_group`` through three
      :meth:`~orpheus.geometry.transformation.RigidMotion.preserves`
      generator checks on the factors (see
      :func:`_derived_product_group`) — measured, not declared;
    * the :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`
      from the construction itself (each polar node is its own level,
      :attr:`~orpheus.numerics.quadrature.rules_sphere.PolarInvariant.SIGNED_MU_Z`);
      :math:`\Sigma` needs nothing stored — it is a query
      (:func:`~orpheus.numerics.symmetry.singular_set`) downstream of
      the nodes.

    Parameters
    ----------
    polar : DiscreteMeasure
        A rule in :math:`\mu = \cos\theta`, nodes shape ``(n_mu,)``
        with :math:`|\mu| \le 1`, carrying an ALGEBRAIC exactness
        claim (e.g.
        :func:`~orpheus.numerics.quadrature.rules_1d.gauss_legendre_on_mu`).
    azimuthal : DiscreteMeasure
        A circle rule on :data:`~orpheus.numerics.measure.SPACE_CIRCLE`,
        nodes shape ``(n_phi, 2)`` of unit points
        :math:`(\cos\varphi, \sin\varphi)`, carrying a TRIGONOMETRIC
        exactness claim (e.g.
        :func:`~orpheus.numerics.quadrature.rules_circle.periodic_trapezoid`).
    atol : float, optional
        Node-matching tolerance for the group derivation's generator
        checks — the same honestly-numerical step as the orbit
        certificate's.

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(n_mu * n_phi, 3)`` on
        ``support=SPACE_SPHERE``, with derived ``invariance_group``
        and ``exactness``.
    LevelStructure
        Per-level indexing for the cylindrical sweep, ordered by the
        stable η-sort (see the tie-break comment at the sort).

    Raises
    ------
    ValueError
        Malformed factors (wrong shapes, non-unit circle nodes, a
        missing claim), mismatched claim systems (via
        :func:`spherical_product_claim`), or a factor pair whose
        derived group family is unrealized (via
        :func:`_derived_product_group`).
    """
    if polar.nodes.ndim != 1:
        raise ValueError(
            f"the polar factor must have scalar μ nodes of shape "
            f"(n_mu,), got shape {polar.nodes.shape} — a spherical "
            f"product's polar factor is a rule in μ = cos(θ)"
        )
    if np.any(np.abs(polar.nodes) > 1.0):
        raise ValueError(
            f"the polar factor's nodes must lie in [-1, 1] "
            f"(μ = cos(θ)); got extremes "
            f"[{polar.nodes.min()}, {polar.nodes.max()}]"
        )
    if azimuthal.support != SPACE_CIRCLE:
        raise ValueError(
            f"the azimuthal factor must be a circle rule on "
            f"{SPACE_CIRCLE!r}, got support {azimuthal.support!r}. An "
            f"interval rule with the same angles is a DIFFERENT object "
            f"carrying a different (algebraic) claim — the confusion "
            f"the exactness carve exists to prevent."
        )
    if azimuthal.nodes.ndim != 2 or azimuthal.nodes.shape[1] != 2:
        raise ValueError(
            f"the azimuthal factor's nodes must be circle POINTS of "
            f"shape (n_phi, 2), got {azimuthal.nodes.shape}"
        )
    radius_defect = np.abs(
        np.linalg.norm(azimuthal.nodes, axis=1) - 1.0
    ).max()
    if radius_defect > atol:
        raise ValueError(
            f"the azimuthal factor's nodes must lie ON the unit "
            f"circle; worst |‖node‖ - 1| = {radius_defect:.3e} — "
            f"off-circle points would produce direction cosines with "
            f"‖Ω‖ ≠ 1"
        )

    # Type-narrowing, not validation: `DiscreteMeasure.exactness` is
    # optional because SOME measures legitimately carry no claim (the
    # direct sum), but a product factor must bring one — the composite's
    # claim is derived from them. Spelled as a raise rather than an
    # `assert` because production asserts vanish under `-O`.
    polar_claim, azimuthal_claim = polar.exactness, azimuthal.exactness
    if polar_claim is None or azimuthal_claim is None:
        raise ValueError(
            f"both factors must carry an exactness claim; got polar="
            f"{polar_claim}, azimuthal={azimuthal_claim}"
        )

    # Derive claim and group BEFORE the O(N) embedding — both refuse on
    # a bad factor pair, and failing fast beats building nodes first.
    claim = spherical_product_claim(polar_claim, azimuthal_claim)
    group = _derived_product_group(polar, azimuthal, atol=atol)

    mu_gl, w_gl = polar.nodes, polar.weights
    cos_phi, sin_phi = azimuthal.nodes[:, 0], azimuthal.nodes[:, 1]
    w_phi = azimuthal.weights
    n_mu = polar.n_points
    n_phi = azimuthal.n_points

    n_total = n_mu * n_phi
    mu_x = np.empty(n_total)
    mu_y = np.empty(n_total)
    mu_z = np.empty(n_total)
    weights = np.empty(n_total)
    level_membership: list[np.ndarray] = []

    idx = 0
    for p in range(n_mu):
        mu_val = mu_gl[p]
        sin_theta = np.sqrt(1.0 - mu_val**2)
        level_idx: list[int] = []
        for m in range(n_phi):
            mu_x[idx] = sin_theta * cos_phi[m]
            mu_y[idx] = sin_theta * sin_phi[m]
            mu_z[idx] = mu_val
            weights[idx] = w_gl[p] * w_phi[m]
            level_idx.append(idx)
            idx += 1
        # MEMBERSHIP ONLY — the ordering belongs to the fiber, and
        # `LevelStructure.from_level_membership` owns it.
        #
        # η = sinθ·cos φ, and `cos φ_m = cos φ_{n_φ−m}` holds BIT-exactly
        # for roots of unity, so a level's η values come in genuine
        # ties: `[M]` only ⌊n_φ/2⌋+1 distinct values among n_φ ordinates
        # (5 of 8 at n_φ=8). Under the `linspace`+`cos` azimuths this
        # rule used until 2026-08-02, round-off manufactured 8 fake
        # distinctions and the sort never saw a tie. This is the
        # 2-to-1-ness the `LevelStructure.level_indices` warning names
        # ("an ordering of the circle modulo the mirror") made literal.
        #
        # Until 2026-08-13 this loop closed with
        # `np.argsort(mu_x[level_arr], kind="stable")`, and the comment
        # here argued that `kind="stable"` was load-bearing. It was — as
        # a REPAIR. Stable keeps the construction order, i.e. increasing
        # φ within an η-tie, which is precisely the fiber key with the φ
        # component left implicit in the loop nesting instead of written
        # down. Naming φ retires the repair: the key is injective, so no
        # tie survives for any algorithm to break. `[M]` the re-route is
        # bit-identical on 8 product configurations (n_μ ∈ {2,3,4,5,6},
        # n_φ ∈ {6,8,10,16,24,32}) and on 5 folded ones — on a FOLDED
        # level η is injective by itself (the T22b theorem the accessor
        # merge rests on), so the fold never saw the tie-break either.
        level_membership.append(np.array(level_idx))

    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=SPACE_SPHERE,
        # COMPUTED from the factors by the three generator checks above
        # — never a declared literal. This module shipped three false
        # symmetry declarations (ERR-072/073/074); a declaration is
        # unfalsifiable, a derivation is refused when its checks fail.
        invariance_group=group,
        # DERIVED from the two factors' own claims, by the theorem that
        # states what a polar × azimuthal embedding into S^2 is exact on.
        exactness=claim,
    )
    structure = LevelStructure.from_level_membership(
        level_membership,
        nodes=nodes,
        level_mu=mu_gl,
        # SIGNED: each polar node is its own level, so the levels run
        # over all of [-1, 1] and a level's fiber is ONE circle. The
        # level-symmetric rule fibers over |mu_z| instead, and its
        # levels carry two.
        polar_invariant=PolarInvariant.SIGNED_MU_Z,
        # The fiber coordinate, as an angle in [0, 2π). Recovered from
        # the circle rule's components rather than kept as a second
        # spelling of the grid — the rule's nodes are POINTS, and the
        # angle is a chart chosen by whoever needs one. The consumers
        # are the arc gates (ω-monotonicity along a folded level) and
        # the curvilinear closure re-pose (arc endpoints, Δω) — and
        # the chart descends to a folded rule by bit-exact SELECTION
        # (`LevelStructure.quotient`), so it is computed exactly once,
        # here.
        #
        # `[M]` The round trip is exact where it matters: against the
        # 2π·m/n_φ grid it reproduces bit-identically at n_φ = 4 and 8
        # (0.0), and to 8.9e-16 at n_φ = 12 — with the ORDER, which is
        # all any consumer reads, preserved at every n_φ tested.
        azimuth=np.tile(
            np.mod(np.arctan2(sin_phi, cos_phi), 2.0 * np.pi), n_mu
        ),
        hemisphere=np.sign(mu_z).astype(np.int64),
    )
    return measure, structure


def product_mu_phi(
    n_mu: int, n_phi: int
) -> tuple[DiscreteMeasure, LevelStructure]:
    r"""Build the GL :math:`\times` equispaced product rule on :math:`S^2`.

    Implements the layout in
    :eq:`product-mu-phi-cosines`. Ordinate ordering matches the
    long-standing SN product-quadrature convention bit-for-bit — the
    one the ``ProductQuadrature`` adapter carried until R-1 Phase A
    detour-C retired the four per-family adapters into classmethod
    factories on the one ``Quadrature`` type: outer loop over the
    :math:`n_\mu`
    Gauss-Legendre nodes (axial cosine), inner loop over the
    :math:`n_\phi` azimuthal samples. Per-level indexing lists are
    sorted by increasing :math:`\eta = \mu_x` to match the
    cylindrical-sweep convention from Bailey-Morel-Chang 2010 Eq. (52).
    ⚠ Eq. (52) states two things; only its ORDER transfers here. Its
    other half - that a cell's radial-cosine measure equals the
    ordinate's weight - is REFUTED on this rule and is not what the
    sort means. See ``docs/theory/methods/sn/curvilinear_one_group.rst`` §bmc-equation-map.

    Weight sum: :math:`\sum_i w_i = 4\pi`.

    Polynomial exactness: the polar factor is exact at
    degree :math:`2 n_\mu - 1` (Stoer-Bulirsch); the azimuthal factor
    is the :math:`n_\phi`-point periodic trapezoid, exact for
    trigonometric polynomials of degree :math:`n_\phi - 1`. We tag
    ``degree_of_exactness = min(2*n_mu - 1, n_phi - 1)`` — the weakest
    of the two, because for general polynomial integrands on
    :math:`S^2` both factors must be exact simultaneously — and the
    bound is **TIGHT**, not conservative: T2/E4 measured the first
    missed degree failing by :math:`1.5\times10^{-2}` while everything
    at or below the claim holds to :math:`4\times10^{-13}`.

    The degree is **derived** from the two factors' own claims by
    :func:`spherical_product_claim`, which carries the theorem. It used
    to be a hand-written ``min()`` over two bare integers; the value is
    unchanged and tight, but it can no longer disagree with the factors
    it is supposed to be about.

    .. note:: What the azimuthal substitution changed (2026-08-02)

       The azimuthal factor was built inline as
       ``np.linspace(0, 2π, n_phi, endpoint=False)`` until this carve. It
       is now
       :func:`~orpheus.numerics.quadrature.rules_circle.periodic_trapezoid`
       at ``shift=NODE_ALIGNED`` — the same grid, named rather than
       spelled, and generated as roots of unity. Three consequences, all
       measured, none cosmetic:

       * :math:`\Sigma` **now has one size.** At even :math:`n_\phi` the
         node set meets the :math:`\xi = 0` axis twice per level, at
         :math:`\varphi = 0` and :math:`\varphi = \pi` — but
         ``np.sin(np.pi)`` is ``1.22e-16``, so under ``linspace`` only
         the first was on it in exact arithmetic. `[M]` At
         ``product_mu_phi(4, 8)``, :math:`|\Sigma|` read **4 by equality
         against 8 by any sane tolerance**; it now reads 8 both ways. A
         fold whose well-posedness condition is :math:`\Sigma =
         \emptyset` cannot be decided on a set whose size depends on how
         you measure it.
       * **The η-ties became real, so the tie-break had to be named.**
         :math:`\cos\varphi_m = \cos\varphi_{n_\phi - m}` holds
         bit-exactly for roots of unity, so a level carries only
         :math:`\lfloor n_\phi/2 \rfloor + 1` distinct :math:`\eta`
         values — `[M]` 5 of 8 at :math:`n_\phi = 8`. Round-off had been
         manufacturing :math:`n_\phi` fake distinctions, which means the
         intra-pair order was previously decided **by noise**. See the
         ``kind="stable"`` comment at the sort.
       * The prose above said "midpoint rule" until 2026-08-02, which was
         false — these are left-endpoints, shift :math:`0`, not shift
         :math:`\tfrac{1}{2}`. The *exactness* sentence beside it was
         right anyway, **because the shift cannot change the degree**.
         What the shift decides is whether a node sits on the mirror
         axis; node-aligned means one does, so :math:`\Sigma \neq
         \emptyset` for every :math:`n_\phi` this rule can be asked for.
         Selecting ``STAGGERED`` instead is the fold's business, and is
         now a substitution rather than a new flag.

    Parameters
    ----------
    n_mu : int
        Number of Gauss-Legendre points in :math:`\mu`. Must be
        :math:`\ge 1`.
    n_phi : int
        Number of equispaced points in :math:`\phi`. Must be
        :math:`\ge 1`.

    Returns
    -------
    DiscreteMeasure
        Nodes shape ``(n_mu * n_phi, 3)``, weights shape
        ``(n_mu * n_phi,)``, on ``support="S^2"``.
        ``invariance_group=SubgroupOfO3.Dnh(n_phi)``,
        ``degree_of_exactness=min(2*n_mu-1, n_phi-1)`` — both DERIVED
        by :func:`spherical_product` from the two factors.
    LevelStructure
        Per-:math:`\mu`-level indexing metadata used by the
        cylindrical SN sweep.

    See Also
    --------
    :meth:`orpheus.numerics.quadrature.Quadrature.product` — the named
    factory SN consumers call. It wraps this measure and holds the
    returned
    :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`
    on its ``level_structure`` field for hot-path access. There is no
    per-family adapter class: the SN-side wrapper this docstring used
    to point at (``orpheus.sn.quadrature.ProductQuadrature``) was
    retired into a classmethod factory on the one ``Quadrature`` type.
    """
    if n_mu < 1:
        raise ValueError(f"product_mu_phi requires n_mu >= 1, got n_mu={n_mu}")
    if n_phi < 1:
        raise ValueError(f"product_mu_phi requires n_phi >= 1, got n_phi={n_phi}")

    # The named family IS a spherical_product of two registered rules:
    # the polar one is the slab rule (one GL source, not a third
    # independent Gauss-Legendre call — consolidated 2026-08-02), and
    # the azimuthal one is the circle rule at NODE_ALIGNED, which
    # reproduces the left-endpoint grid this rule has always used. The
    # nodes are roots of unity, so cos/sin are exact on the axes — what
    # makes Sigma = {xi = 0} a set of a definite size (see the note in
    # the docstring above).
    return spherical_product(
        gauss_legendre_on_mu(n_mu),
        periodic_trapezoid(n_phi, shift=NODE_ALIGNED),
    )


__all__ = ["product_mu_phi", "spherical_product", "spherical_product_claim"]

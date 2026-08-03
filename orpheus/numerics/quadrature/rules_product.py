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
Together these generate :math:`D_{n_\phi h}` — and
:func:`~orpheus.numerics.symmetry.maximal_invariance_groups`
**computes** exactly that from the nodes, so the tag above is pinned
by a permanent test rather than merely asserted. It is not verified
at construction: the computed check runs in ``tests/``, not here.

.. caution::

   The agreement holds to ~1e-16, not bit-exactly, and the match window
   (``symmetry._NODE_WINDOW_FACTOR``) is what absorbs the difference.
   Lebedev's :math:`O_h` claim, by contrast, is exact on both sides
   because signed permutations are exact in IEEE.

   Issue #325 named **two** halves for closing this gap, and as of
   2026-08-02 the **generator half is done**: the :math:`\phi` grid is
   :mod:`~orpheus.numerics.roots_of_unity`-generated, so this side of
   the check no longer evaluates ``cos``/``sin`` at a sampled angle.
   The remaining half is the **checker's** own :math:`C_n` and
   :math:`\sigma_v` operators, which still do — so the window is still
   load-bearing, now for one reason instead of two.

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

* Bailey, T.S., Adams, M.L., Yang, B., and Zika, M.R. (2009). "A
  piecewise linear discontinuous finite element spatial discretization
  of the transport equation." *Annals of Nuclear Energy* **35**,
  1929-1936. Eq. 50 — α-recursion convention used by the cylindrical
  sweep that consumes this product rule.
"""

from __future__ import annotations

import numpy as np

from ..exactness import (
    UNIFORM_ON_SPHERE,
    ExactnessClaim,
    OrthogonalSystem,
)
from ..measure import SPACE_SPHERE, DiscreteMeasure
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


def product_mu_phi(
    n_mu: int, n_phi: int
) -> tuple[DiscreteMeasure, LevelStructure]:
    r"""Build the GL :math:`\times` equispaced product rule on :math:`S^2`.

    Implements the layout in
    :eq:`product-mu-phi-cosines`. Ordinate ordering matches the
    long-standing :class:`orpheus.sn.quadrature.ProductQuadrature`
    convention bit-for-bit: outer loop over the :math:`n_\mu`
    Gauss-Legendre nodes (axial cosine), inner loop over the
    :math:`n_\phi` azimuthal samples. Per-level indexing lists are
    sorted by increasing :math:`\eta = \mu_x` to match the
    cylindrical-sweep convention from Bailey et al. (2009) Eq. 50.

    Weight sum: :math:`\sum_i w_i = 4\pi`.

    Polynomial exactness: the polar factor is exact at
    degree :math:`2 n_\mu - 1` (Stoer-Bulirsch); the azimuthal factor
    is the :math:`n_\phi`-point periodic trapezoid, exact for
    trigonometric polynomials of degree :math:`n_\phi - 1`. We tag the
    conservative ``degree_of_exactness = min(2*n_mu - 1, n_phi - 1)``
    as the weakest of the two — for general polynomial integrands on
    :math:`S^2`, both factors must be exact simultaneously.

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
        ``(n_mu * n_phi,)``, on ``space="S^2"``.
        ``invariance_group=SubgroupOfO3.Dnh(n_phi)``,
        ``degree_of_exactness=min(2*n_mu-1, n_phi-1)``.
    LevelStructure
        Per-:math:`\mu`-level indexing metadata used by the
        cylindrical SN sweep.

    See Also
    --------
    :class:`orpheus.sn.quadrature.ProductQuadrature` — the SN-side
    adapter caching this rule's outputs for hot-path access.
    """
    if n_mu < 1:
        raise ValueError(f"product_mu_phi requires n_mu >= 1, got n_mu={n_mu}")
    if n_phi < 1:
        raise ValueError(f"product_mu_phi requires n_phi >= 1, got n_phi={n_phi}")

    # GL points in μ = cos(θ). The polar factor IS the slab rule — same
    # construction, one source — rather than a third independent call to
    # a Gauss-Legendre routine. Consolidated 2026-08-02; it had been
    # `np.polynomial.legendre.leggauss(n_mu)`, which the audit found as
    # a surviving twin of the two rules the campaign was retiring.
    polar = gauss_legendre_on_mu(n_mu)
    mu_gl, w_gl = polar.nodes, polar.weights

    # The azimuthal factor IS the registered circle rule — same reason as
    # the polar one above. NODE_ALIGNED reproduces the left-endpoint grid
    # this rule has always used; it is now a named shift rather than a
    # `linspace` literal, which is what lets a fold SELECT the staggered
    # one instead of a flag existing.
    #
    # The nodes are roots of unity, so cos/sin are exact on the axes.
    # That is not cosmetic here: it is what makes Sigma = {xi = 0} a set
    # of a definite size (see the note in the docstring above).
    azimuthal = periodic_trapezoid(n_phi, shift=NODE_ALIGNED)
    cos_phi, sin_phi = azimuthal.nodes[:, 0], azimuthal.nodes[:, 1]
    w_phi = azimuthal.weights

    # Type-narrowing, not validation: `DiscreteMeasure.exactness` is
    # optional because SOME measures legitimately carry no claim (the
    # direct sum), but both factors here are built two lines up by rules
    # whose return contract includes one. Spelled as a raise rather than
    # an `assert` because production asserts vanish under `-O`.
    polar_claim, azimuthal_claim = polar.exactness, azimuthal.exactness
    if polar_claim is None or azimuthal_claim is None:  # pragma: no cover
        raise ValueError(
            f"both factors must carry an exactness claim; got polar="
            f"{polar_claim}, azimuthal={azimuthal_claim}"
        )

    n_total = n_mu * n_phi
    mu_x = np.empty(n_total)
    mu_y = np.empty(n_total)
    mu_z = np.empty(n_total)
    weights = np.empty(n_total)
    level_indices: list[np.ndarray] = []

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
        # Sort by increasing η (mu_x) for cylindrical sweep convention.
        #
        # `kind="stable"` is load-bearing, and became so the moment the
        # azimuths turned exact. η = sinθ·cos φ, and `cos φ_m = cos
        # φ_{n_φ−m}` holds BIT-exactly for roots of unity — so a level's
        # η values come in genuine ties, `[M]` only ⌊n_φ/2⌋+1 distinct
        # values among n_φ ordinates (5 of 8 at n_φ=8). Under the
        # `linspace`+`cos` azimuths this rule used until 2026-08-02,
        # round-off manufactured 8 fake distinctions and the sort never
        # saw a tie.
        #
        # This is the 2-to-1-ness `LevelStructure.fiber` already names
        # ("ordered by a projection that is 2-to-1 on it") made literal.
        # With real ties the tie-break must be NAMED rather than left to
        # the sort algorithm: stable keeps the construction order, i.e.
        # increasing φ within an η-tie.
        level_arr = np.array(level_idx)
        order = np.argsort(mu_x[level_arr], kind="stable")
        level_indices.append(level_arr[order])

    nodes = np.column_stack([mu_x, mu_y, mu_z])  # (N, 3)
    measure = DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=SPACE_SPHERE,
        # D_{n_φ h}, NOT SO(2). The φ grid is a FINITE set of n_φ
        # equispaced azimuths, so no continuous rotation preserves it —
        # only the n_φ-fold cyclic rotations, the vertical mirrors at
        # kπ/n_φ, and (because the GL μ nodes are symmetric) σ_h. The
        # rule carried ``SO2`` until 2026-08-02: a claim no finite point
        # set on S² can satisfy. See ``maximal_invariance_groups``, which
        # computes this group from the nodes and pins the declaration.
        invariance_group=SubgroupOfO3.Dnh(n_phi),
        # DERIVED from the two factors' own claims, by the theorem that
        # states what a polar × azimuthal embedding into S^2 is exact on.
        # The value is unchanged — `min(2n_μ−1, n_φ−1)`, tight — but it is
        # no longer a hand-written `min()` over two bare integers that
        # happen to match the factors' degrees. It now cannot drift from
        # them, and the systems are checked rather than assumed.
        exactness=spherical_product_claim(polar_claim, azimuthal_claim),
    )
    structure = LevelStructure(
        n_levels=n_mu,
        level_indices=level_indices,
        level_mu=mu_gl,
        # SIGNED: each GL polar node is its own level, so the levels run
        # over all of [-1, 1] and a level's fiber is ONE circle. The
        # level-symmetric rule fibers over |mu_z| instead, and its
        # levels carry two.
        polar_invariant=PolarInvariant.SIGNED_MU_Z,
        # The fiber coordinate, as an angle in [0, 2π). Recovered from
        # the circle rule's components rather than kept as a second
        # spelling of the grid — the rule's nodes are POINTS, and the
        # angle is a chart chosen by whoever needs one. `fiber()` needs
        # one, because it orders by it.
        #
        # `[M]` The round trip is exact where it matters: against the
        # 2π·m/n_φ grid it reproduces bit-identically at n_φ = 4 and 8
        # (0.0), and to 8.9e-16 at n_φ = 12 — with the ORDER, which is
        # all `fiber()` consumes, preserved at every n_φ tested.
        azimuth=np.tile(
            np.mod(np.arctan2(sin_phi, cos_phi), 2.0 * np.pi), n_mu
        ),
        hemisphere=np.sign(mu_z).astype(np.int64),
    )
    return measure, structure


__all__ = ["product_mu_phi", "spherical_product_claim"]

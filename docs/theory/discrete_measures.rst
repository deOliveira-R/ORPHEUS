.. _discrete-measures:

============================================
Discrete Measures and Quadrature Composition
============================================

.. note::

   **Stub page** — full narrative pending under Issue 18 of
   ``.claude/plans/sn_reshape.md``. This stub installs the
   load-bearing equation labels (:eq:`discrete-measure-integrate`,
   :eq:`discrete-measure-pushforward`,
   :eq:`bundle-measure-disintegration`) and the cross-references
   that downstream theory pages and tests rely on.

Key Facts
=========

- A **quadrature rule** is a discrete measure
  :math:`\mu = \sum_i w_i \, \delta_{x_i}` on a measurable space
  :math:`(\mathcal{X}, \Sigma)`. Lebesgue integration of a test
  function against :math:`\mu` is the familiar quadrature sum
  (:eq:`discrete-measure-integrate`).
- Four canonical composition operations are exposed by the project's
  primitive :class:`~orpheus.numerics.measure.DiscreteMeasure`:

  - **Tensor product** ``μ * ν`` — product measure
    :math:`\mu \otimes \nu` on
    :math:`\mathcal{X} \times \mathcal{Y}`.
  - **Direct sum** ``μ + ν`` — concatenation of atoms on a shared
    space; the foundation of composite (panelled) quadrature.
  - **Pushforward** ``μ.pushforward(φ)`` — image measure
    :math:`\varphi_* \mu` (:eq:`discrete-measure-pushforward`).
  - **Restriction** ``μ.restrict(E)`` — indicator-multiplication
    :math:`\mathbf{1}_E \cdot \mu`; supports half-range SN sweeps
    and bundle cuts.

- :class:`~orpheus.numerics.measure.BundleMeasure` carries the
  **disintegration** :eq:`bundle-measure-disintegration`: a base
  measure on :math:`\mathcal{B}` plus, for each base point, a fiber
  measure on :math:`\pi^{-1}(b)`. SN does not consume bundles in the
  Wave A campaign; MoC ray-bundle quadratures (Wave 2) will.
- The space tag is a runtime ``str`` — Python generics over a
  measurable-space type are not expressive enough without runtime
  overhead, so the project uses opaque string tags for sanity checks
  and documentation only.


Definitions
===========

A discrete (atomic) measure on a measurable space
:math:`(\mathcal{X}, \Sigma)` is a finite sum of weighted Dirac
deltas:

.. math::
   :label: discrete-measure-definition

   \mu = \sum_{i=1}^{N} w_i \, \delta_{x_i},
   \qquad x_i \in \mathcal{X},\quad w_i \in \mathbb{R}.

.. vv-status: discrete-measure-definition documented

The Lebesgue integral of an integrable test function :math:`f`
against :math:`\mu` collapses to a finite sum,

.. math::
   :label: discrete-measure-integrate

   \int_{\mathcal{X}} f \, d\mu \;=\; \sum_{i=1}^{N} w_i \, f(x_i),

which is the **quadrature** identity. Polynomial-exact rules
(Gauss-Legendre, Gauss-Chebyshev, …) are characterised by the
maximum degree :math:`p` such that
:eq:`discrete-measure-integrate` is exact for every polynomial of
degree :math:`\le p` — the field's
:attr:`~orpheus.numerics.measure.DiscreteMeasure.degree_of_exactness`
records this.

For a measurable map
:math:`\varphi : \mathcal{X} \to \mathcal{Y}`, the **pushforward**
or image measure :math:`\varphi_* \mu` satisfies the
change-of-variables identity

.. math::
   :label: discrete-measure-pushforward

   \int_{\mathcal{Y}} g \, d(\varphi_* \mu)
   \;=\; \int_{\mathcal{X}} (g \circ \varphi) \, d\mu
   \;=\; \sum_i w_i \, g(\varphi(x_i)).

.. (vv-status rationale) discrete-measure-pushforward: Verified
   transitively by the foundation test
   ``test_pushforward_invertible_map_change_of_variables`` —
   integrates :math:`x^2` against the φ-image of an equispaced
   measure on :math:`[0, 1]` under :math:`\varphi(x) = x^2` and
   confirms convergence to :math:`\int_0^1 x^4 \, dx = 1/5`,
   exactly the change-of-variables identity in this label. Tagged
   ``foundation`` rather than carrying ``verifies(...)`` because
   the identity is a software invariant of the discrete pushforward
   operation, not a physics-equation claim with an L0..L3 ladder
   slot.
.. vv-status: discrete-measure-pushforward documented

The pushforward in
:class:`~orpheus.numerics.measure.DiscreteMeasure` is
**φ-image semantics**: weights are preserved verbatim. If the
caller wants Radon-Nikodym semantics
:math:`g(y)\,dy = h(x)\,dx`, the Jacobian
:math:`|\det D\varphi|^{-1}` must be applied to the weights
explicitly.

For a fibered space :math:`\pi : \mathcal{X} \to \mathcal{B}`, the
**disintegration theorem** (Bourbaki 1969, Intégration VI §3) gives

.. math::
   :label: bundle-measure-disintegration

   \int_\mathcal{X} f \, d\mu
   \;=\; \int_\mathcal{B} \!\left[ \int_{\pi^{-1}(b)} f \, d\mu_b
   \right] d\mu_\mathcal{B}(b),

.. (vv-status rationale) bundle-measure-disintegration: Verified
   transitively by the foundation tests
   ``test_bundle_measure_separable_smoke`` (independent fiber, base
   :math:`\otimes` fiber agrees with explicit nested sum) and
   ``test_bundle_measure_varying_fiber`` (fiber depends on base
   point, dispatch path is exercised). MoC's ray-bundle quadratures
   in Wave 2 of the SN-reshape campaign will provide the physics-grade
   verification once the abstraction is consumed.
.. vv-status: bundle-measure-disintegration documented

implemented in
:class:`~orpheus.numerics.measure.BundleMeasure` as a base measure
plus a lazy fiber-measure factory.


Composition algebra — metadata propagation
==========================================

The four composition operations propagate the
:class:`~orpheus.numerics.measure.DiscreteMeasure` metadata as
follows. Where a field is marked **dropped**, the operation
generally cannot guarantee the field is still meaningful, and the
result carries ``None`` until the caller re-asserts it via
:meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`.

.. list-table::
   :header-rows: 1
   :widths: 18 18 16 16 16 16

   * - Operation
     - ``space``
     - ``n_points``
     - ``dim``
     - ``degree_of_exactness``
     - ``invariance_group``
   * - tensor product ``μ * ν``
     - ``f"{μ.space} × {ν.space}"``
     - :math:`N_\mu \cdot N_\nu`
     - :math:`d_\mu + d_\nu`
     - :math:`\min(p_\mu, p_\nu)` if both set, else dropped
     - dropped (Issue 3 may refine when both factors carry
       compatible groups)
   * - direct sum ``μ + ν``
     - shared (raises if mismatched)
     - :math:`N_\mu + N_\nu`
     - shared (raises if mismatched)
     - :math:`\min(p_\mu, p_\nu)` if both set, else dropped
     - dropped (concatenation generally breaks invariance)
   * - pushforward ``μ.pushforward(φ)``
     - explicit ``new_space=`` argument or ``f"φ_*({μ.space})"``
     - unchanged
     - inferred from ``φ``'s output shape
     - dropped (φ-image of a polynomial-exact rule is not
       polynomial-exact in general — undefined unless ``φ`` is
       affine)
     - dropped (φ generally does not preserve a group action)
   * - restriction ``μ.restrict(E)``
     - unchanged
     - ``mask.sum()``
     - unchanged
     - dropped (restricted rule is not Gauss on the cut domain)
     - dropped (caller may re-tag if :math:`E` is invariant)

Be honest about what is dropped: the field becoming ``None`` is the
correct behaviour, not a degradation. If a downstream consumer
needs the field, it must re-establish the claim and stamp it back
with :meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`.


1-D primitive constructors
==========================

The module ships three 1-D primitives:

- :func:`~orpheus.numerics.measure.gauss_legendre`
  ``(n)`` — Gauss-Legendre on :math:`[-1, 1]`,
  ``degree_of_exactness = 2n - 1``, weights sum to 2.
- :func:`~orpheus.numerics.measure.gauss_chebyshev`
  ``(n)`` — Gauss-Chebyshev (first kind) on :math:`[-1, 1]` with
  weight :math:`(1 - x^2)^{-1/2}`, ``degree_of_exactness = 2n - 1``
  (in the weighted sense), weights sum to :math:`\pi`.
- :func:`~orpheus.numerics.measure.equispaced`
  ``(a, b, n)`` — midpoint rule on :math:`[a, b]`,
  ``degree_of_exactness = 1``, weights sum to :math:`b - a`.

Higher-dimensional rules are built by composition:
``μ_S2 = gauss_legendre(n_mu) * equispaced(0.0, 2*np.pi, n_phi)``
gives a product quadrature on
:math:`\mu \in [-1, 1] \times \varphi \in [0, 2\pi)` —
the standard polar-azimuthal split for transport on the unit
sphere.


Forward references
==================

This page is consumed progressively. The current stub installs the
labels and the composition algebra; the following pieces will land
in subsequent issues of the SN reshape campaign:

- **Issue 3** (``module:numerics``) will populate the
  ``invariance_group`` field with a tag for the symmetry group
  under which a measure is invariant — :math:`O(3)`, finite
  subgroups, the symmetric group, etc. — and document the
  composition rule for invariance under tensor product. The
  ``invariance_group`` field is reserved here for that integration
  point.
- **Issue 4** (``module:sn``) will refactor
  :class:`~orpheus.sn.quadrature.AngularQuadrature` so that
  :class:`~orpheus.sn.quadrature.GaussLegendre1D`,
  :class:`~orpheus.sn.quadrature.LebedevSphere`,
  :class:`~orpheus.sn.quadrature.LevelSymmetricSN`, and
  :class:`~orpheus.sn.quadrature.ProductQuadrature` are thin
  adapters wrapping a
  :class:`~orpheus.numerics.measure.DiscreteMeasure` rather than
  carrying their own ``mu_x``/``mu_y``/``weights`` fields. This
  hides the composition behind the existing SN API while exposing
  it to consumers (such as the eigenvalue solver tests, the MOC
  port, and Wave-3 sensitivity sweeps) that now want to reason
  about quadratures as measures.
- **Wave 2 (MoC migration)** will consume
  :class:`~orpheus.numerics.measure.BundleMeasure` for ray-bundle
  quadratures whose fibers (the parallel rays through the geometry)
  vary with the base direction. The bundle abstraction is shipped
  in Phase 0 because retrofitting it later would force a rewrite
  of MoC's quadrature handling.

The narrative in :ref:`theory-reference-solvers` and
:ref:`theory-discrete-ordinates` will progressively replace any
quadrature-specific wording with cross-references to
:eq:`discrete-measure-integrate` and the composition algebra
above.


References
==========

* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. Springer. §3.6 (Gauss quadrature);
  Theorem 3.6.20 (polynomial exactness).
* Trefethen, L.N. (2008). "Is Gauss quadrature better than
  Clenshaw-Curtis?" *SIAM Review* **50**, 67-87. The measure-theoretic
  perspective on spectral integration rules.
* Bourbaki, N. (1969). *Éléments de mathématique. Intégration*,
  Chapter VI §3 (disintegration of measures).

.. _theory-glossary:

Glossary
========

Canonical reactor-physics and transport terms used across the theory pages.
Cross-reference a term from any page with the ``:term:`` role (e.g.
``:term:`optical thickness```) — this both links the reader to the definition
and creates a page → concept edge in the Nexus knowledge graph, populating the
``TERM`` node layer for concept-level retrieval (:ref:`#231 <theory-glossary>`).

This glossary is seeded from the SN prototype (issue #231) and grows as each
theory page adopts ``:term:`` references. Definitions double as concept-level
embedding anchors, so each stands alone in query vocabulary.

.. glossary::
   :sorted:

   albedo
      The fraction of neutrons incident on a surface that are returned to the
      domain: :math:`\alpha = J^-/J^+`, the ratio of the returning to the
      outgoing partial current. :math:`\alpha = 0` is a :term:`vacuum boundary
      condition`; :math:`\alpha = 1` is a perfect reflector.

   angular flux
      :math:`\psi(\mathbf{r}, \mathbf{\Omega}, E)` — the expected track length
      per unit volume, solid angle, energy, and time of neutrons at position
      :math:`\mathbf{r}` travelling in direction :math:`\mathbf{\Omega}` with
      energy :math:`E`. The fundamental unknown of the transport equation.

   scalar flux
      :math:`\phi(\mathbf{r}, E) = \int \psi(\mathbf{r}, \mathbf{\Omega}, E)\,
      d\mathbf{\Omega}` — the :term:`angular flux` integrated over all
      directions (the zeroth angular moment). Reaction rates are
      :math:`\phi\,\Sigma`.

   diamond difference
      The spatial closure :math:`\bar{\psi} = \tfrac{1}{2}(\psi_{\text{in}} +
      \psi_{\text{out}})` relating a cell-average :term:`angular flux` to its
      incoming and outgoing face values. Second-order accurate but not positive;
      the :math:`\tau = \tfrac{1}{2}` case of the :term:`weighted diamond
      difference`.

   weighted diamond difference
      The generalization :math:`\bar{\psi} = \tau\,\psi_{\text{out}} +
      (1-\tau)\,\psi_{\text{in}}` with a weight :math:`\tau \in [0, 1]`;
      :math:`\tau = \tfrac{1}{2}` recovers :term:`diamond difference`. In
      curvilinear geometry the Morel--Montry weight sets :math:`\tau` per
      :term:`ordinate` to stay exact for a flux linear in :math:`\mu`.

   lethargy
      :math:`u = \ln(E_0/E)` — a dimensionless logarithmic energy variable that
      increases as energy decreases. Group widths in lethargy,
      :math:`\Delta u = \ln(E_g/E_{g+1})`, are the natural measure for the
      slowing-down density.

   scattering ratio
      :math:`c = \Sigma_s/\Sigma_t` — the mean number of secondaries per
      collision (equivalently, the probability that a collision is a scatter
      rather than an absorption). :math:`c = 1` is a purely scattering
      (conservative) medium and sets the source-iteration spectral radius
      :math:`\rho = c`.

   optical thickness
      :math:`\tau = \int \Sigma_t\, ds` — path length measured in mean free
      paths; the dimensionless argument of the exponential attenuation
      :math:`e^{-\tau}`.

   ordinate
      A discrete direction :math:`\mathbf{\Omega}_n`, carrying a
      :term:`quadrature` weight :math:`w_n`, at which the :term:`angular flux`
      is evaluated in the discrete-ordinates (:math:`S_N`) method.

   quadrature
      A set of :term:`ordinate`\ s and weights :math:`\{\mathbf{\Omega}_n,
      w_n\}` approximating the angular integral :math:`\int f\, d\mathbf{\Omega}
      \approx \sum_n w_n\, f(\mathbf{\Omega}_n)`. The weights sum to the total
      solid angle (:math:`4\pi` in three dimensions).

   sweep
      The ordered, direction-by-direction spatial solve that inverts the
      streaming--collision operator :math:`(L + C)^{-1}` for a fixed source,
      marching cell-by-cell from the inflow boundary along each
      :term:`ordinate`.

   starting direction
      In curvilinear (spherical / cylindrical) :math:`S_N`, the tangent
      :term:`ordinate` :math:`\mu = \pm 1` at which the angular-redistribution
      recurrence is seeded. The Carlson starting-direction sweep supplies the
      :math:`\psi_{1/2}` half-angle flux the redistribution needs.

   reflective boundary condition
      A specular boundary that returns each incident :term:`ordinate` into its
      mirror-image direction across the surface normal; conserves neutrons
      (:term:`albedo` 1).

   vacuum boundary condition
      A boundary that returns nothing: the incoming :term:`angular flux` is zero
      for all inward directions (:term:`albedo` 0). Models a non-reentrant
      surface.

   white boundary condition
      An isotropic-return boundary: incident neutrons are re-emitted with a
      direction-independent distribution whose magnitude conserves the incident
      current. Used for approximate reflective symmetry on curved surfaces.

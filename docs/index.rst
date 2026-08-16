ORPHEUS — Open Reactor Physics Educational University System
=============================================================

Neutron-transport theory and verified solvers, documented as one
corpus: full derivations, design rationale, conventions, gotchas, and
what verifies each equation.

**Where to start.** :doc:`theory/index` is the canonical entry — its
part table routes by task (importing an equation, touching a solver,
writing a reference). The S\ :sub:`N` sub-book is the most developed
method treatment (:doc:`theory/methods/sn/index`, with reading tracks
and a symptom → chapter diagnostic table); the verification part —
principles, harness, and the auto-generated V&V matrix — lives at
:doc:`theory/verification/index`.

.. toctree::
   :maxdepth: 2
   :caption: Theory & Derivations

   theory/index

.. toctree::
   :maxdepth: 2
   :caption: Architecture

   architecture/index

.. toctree::
   :maxdepth: 2
   :caption: Development

   development

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/numerics
   api/transport
   api/data
   api/geometry
   api/homogeneous
   api/discrete_ordinates
   api/collision_probability
   api/method_of_characteristics
   api/monte_carlo
   api/diffusion_1d
   api/fuel_behaviour
   api/thermal_hydraulics
   api/reactor_kinetics
   api/derivations


Knowledge Graph
================

.. The href below is a raw path into the build output, so nothing verifies
   it — Sphinx does not check relative hyperlinks, at any severity. It must
   mirror ``[graph].output`` in ``.nexus/config.toml``; if that setting
   moves, this link moves with it or it 404s silently.

`Open interactive graph explorer <graph/graph.html>`_ — visualize the full code + documentation connectivity graph with filtering, search, and node inspection.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Specification (``specification``)
==================================

The :mod:`orpheus.specification` package holds the reference
specification: the question a reference answers, written down with the
materials and, unless the problem is the infinite medium, the geometry it
is asked of. It is the key a reference cache stores an answer under, so
it is a content value admitted in a canonical form at construction. The
layer of the posing filtration the question is posed at is the type:

* :class:`~orpheus.specification.specification.InfiniteMediumSpecification`
  ``(material_id, mixture, question)`` is posed on energy alone, the
  infinite medium (:ref:`infinite-medium-definition`), and has no
  geometry field;
* :class:`~orpheus.specification.specification.GeometrySpecification`
  ``(materials, geometry, question)`` is posed on a finite
  :class:`~orpheus.geometry.structured_geometry.StructuredGeometry` with
  its boundary laws, and keeps only the materials the geometry assigns.

``Specification`` is the closed union of the two, and ``Coordinate =
CellCoefficient | GeometryExtent`` the closed set of keys a question is
resolved to (:mod:`orpheus.data.cells`, :mod:`orpheus.geometry.extent`).

Whether an observable (:mod:`orpheus.numerics.observable`) can be posed on
a specification is decided once, by
:func:`~orpheus.specification.specification.admit_observable`: a flux
integral's weight must fit the problem as a question's datum must, a
ratio's two operands must each fit, a point value's group must be one of
the problem's and its position must lie on the geometry (so the infinite
medium refuses it), and an eigenvalue exists only for an eigen question.
A published solution calls it on every observable it prints
(:mod:`orpheus.reference.published`).

The package is input-tier: it imports :mod:`orpheus.data`,
:mod:`orpheus.geometry` and :mod:`orpheus.numerics` and nothing above
them, and :mod:`orpheus.derivations` may import it
(:ref:`architecture-layering`). The theory, with why the layer is the
type, the canonical form, the coordinates, the datum's fit, what is
deferred and the gates, is :ref:`structured-geometry-specification`.

.. automodule:: orpheus.specification

.. automodule:: orpheus.specification.specification

.. autoclass:: orpheus.specification.specification.InfiniteMediumSpecification
   :members: materials, n_groups, n_regions, unreadable, resolve_extent
   :undoc-members:

.. autoclass:: orpheus.specification.specification.GeometrySpecification
   :members: n_groups, n_regions, unreadable, resolve_extent
   :undoc-members:

.. autodata:: orpheus.specification.specification.Specification

.. autodata:: orpheus.specification.specification.Coordinate

.. autofunction:: orpheus.specification.specification.admit_observable

Values (``values``)
===================

The :mod:`orpheus.values` package is the bottom of the package graph: what a
value is, and when two values are the same. A value type is a set with an
equality, and the package holds the two halves. Every layer that admits a
value imports them from here, and the package imports nothing from
``orpheus`` (the gate is ``tests/gates/test_layer_imports.py``; the reason for
the package and its name is ruling 6 of the layer DAG campaign, and the layer
it occupies is on :doc:`/architecture/layering`). Both modules were
``orpheus.numerics.scalars`` and ``orpheus.numerics.content`` until
2026-10-22; content identity does not key on a class's module (ruling 13), so
the move changed no digest.

.. automodule:: orpheus.values

The set
-------

Scalars — :mod:`orpheus.values.scalars`, the one definition of "a real
number" that the content encoder, the data layer, the geometry, the mesh
and the numerics values all parse through: NaN refused, ``-0.0`` made
``+0.0``, an integer beyond :math:`2^{53}` and an overflow refused with
their key. It was ``orpheus.geometry.scalars`` until #559. The theory,
with what each function adds and why the encoder admits ``bool`` while
the parsers refuse it, is :ref:`structured-geometry-one-real-parser`.

.. automodule:: orpheus.values.scalars
   :members:

The equality
------------

Content identity — :mod:`orpheus.values.content`, the one encoder
behind every persistent key and every content ``==`` and ``hash``. Its
module docstring is the design record; the theory, with the canonical
form of each kind of value and its reason, the schema tag, the two
refusals and the gates, is :ref:`structured-geometry-content-identity`.

.. py:module:: orpheus.values.content

.. autofunction:: orpheus.values.content.encode

.. autofunction:: orpheus.values.content.encode_exact

.. autofunction:: orpheus.values.content.content_digest

.. autofunction:: orpheus.values.content.name_digest

.. autoexception:: orpheus.values.content.ContentlessError

.. autoclass:: orpheus.values.content.ContentIdentity
   :members: content_parts, content_digest

.. autoclass:: orpheus.values.content.FrozenMapping

.. autofunction:: orpheus.values.content.constructor_arguments

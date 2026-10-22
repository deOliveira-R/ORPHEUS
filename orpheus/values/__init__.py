r"""Values: what a value is, and when two values are the same.

A value type is a set with an equality. This package holds the two halves,
and every layer of ORPHEUS that admits a value imports them from here:

* :mod:`orpheus.values.scalars` is the set: the one definition of a real
  number, a finite real, an integer and a positive quantity (NaN refused,
  ``-0.0`` made ``+0.0``, an integer a double would round refused);
* :mod:`orpheus.values.content` is the equality: two values are the same
  when their content is the same, one encoder and one digest, from which
  :class:`~orpheus.values.content.ContentIdentity` derives ``==`` and
  ``hash``.

It is the lowest layer of the package graph and imports nothing from
``orpheus`` (``tests/gates/test_layer_imports.py``). A module joins only if
it defines what a value is or when two values are equal; the reason for the
package, and for its name, is ruling 6 of ``.claude/plans/layer_dag.md``.
"""

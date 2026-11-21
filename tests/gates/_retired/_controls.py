r"""Planted sources that spell REAL retired names, for the mutations of sub-unit ``angular-mass`` (``ontology_w0_spec.md`` 5.7.3 E6 and 5.7.5 G3).

This file lives in ``tests/gates/_retired/`` because that package is excluded from the retirement engine's scan: a planted
source spelling ``CONTAINS``, ``__add__`` or ``__radd__`` of ``AngularSourceSink`` and ``BulkField`` in an ordinary test file would redden
its own rows (the string leg parses every string constant of a scanned file). Nothing here is imported by a gate that runs today; the
reviewers' copy-based mutations read these constants.

* ``ANGULAR_SOURCE_SINK_ARM``: the cross-class arm of ``orpheus/transport/fields/angular_source_sink.py`` at the base tip
  ``954c4a4f``, verbatim, from the section comment ``# -- Algebra extensions (over Field)`` through ``return NotImplemented`` of
  ``__radd__``. The mutant also restores three imports (``replace``, ``BulkField``, ``ScalarSourceSink``), or it fails at import
  with ``NameError`` and is never red evidence.
* ``BULK_FIELD_CONTAINS``: ``BulkField.CONTAINS`` and ``BulkField.__add__`` of ``orpheus/transport/fields/_bases.py``, verbatim.
* ``PLANTED_BULK_FIELD_CONTAINS``: the one-line plant of the engine mutation (insert after the end of the ``BulkField`` class
  docstring; the ``BulkField.CONTAINS`` row must be in the red set of the walk on that root).

Each constant is the base text of the file it was cut from (it occurs exactly once there, measured at the base tip by slicing the two
files and comparing the slices byte for byte); after the unit lands the text is gone from the tree, which is the point of a plant.
"""

ANGULAR_SOURCE_SINK_ARM = r'''    # ── Algebra extensions (over Field) ──────────────────────────────

    #: The isotropic source is the subspace where every ordinate carries the
    #: same value (``iso -> 1 (x) iso``), so this class answers ``iso + aniso``.
    CONTAINS: ClassVar[tuple[type[BulkField], ...]] = (ScalarSourceSink,)

    def __add__(self, other):
        r"""Add a partner source, dispatched on the partner's type.

        * An :class:`AngularSourceSink` partner: the within-class add, the
          algebra :class:`~orpheus.numerics.field.Field` provides.
        * A :class:`ScalarSourceSink` partner: the canonical subspace
          containment. The iso operand is broadcast across the Ω axis (the
          structural map ``iso → 1 ⊗ iso``) and added into this per-ordinate
          operand's space, and the result is an :class:`AngularSourceSink`.
          The sum is commutative: ``iso + aniso`` reaches :meth:`__radd__`,
          because ``BulkField.__add__`` defers to the containing class.

        Coherence is the space relation: the iso operand's space must be this
        operand's non-angular marginal (drop axis 0), compared axis-wise by
        content.

        The module docstring gives the justification (refined Issue #207:
        cross-class dunders are permitted through a canonical subspace
        containment). The dunder is deliberately untyped: its cross-class arm
        returns the larger space's class, which no static spelling reconciles
        with :meth:`Field.__add__`'s ``(T, T) -> T`` contract (#288).
        """
        if isinstance(other, ScalarSourceSink):
            if other.space.axes is None or self.space.axes is None or self.space.axes[1:] != other.space.axes:
                raise ValueError(
                    "ScalarSourceSink + AngularSourceSink requires the iso "
                    "operand's space to be the angular operand's non-angular "
                    "marginal (energy ⊗ spatial[⊗ moment] content match); "
                    f"got {other.space!r} vs the marginal of {self.space!r}."
                )
            return replace(self, values=other.values[None] + self.values)
        return super().__add__(other)

    def __radd__(self, other):
        r"""``iso + aniso``: the containment sum, commutative (:meth:`__add__`)."""
        if isinstance(other, ScalarSourceSink):
            return self.__add__(other)
        return NotImplemented
'''

BULK_FIELD_CONTAINS = r'''    # ── Algebra extension (over Field) ───────────────────────────────

    #: The bulk-field classes whose space embeds canonically into this class's
    #: (an isotropic source is the angular source constant across the
    #: ordinates: ``iso -> 1 (x) iso``). A sum with one of them is answered by
    #: THIS class, the containing one, so the embedding is written once, in the
    #: class of the larger space, and the smaller class names no larger one.
    CONTAINS: ClassVar[tuple[type["BulkField"], ...]] = ()

    def __add__(self, other):
        """Defer to a partner whose space contains this one's (its reflected add answers); else the base algebra.

        Deliberately untyped, like the containing class's dunder: the
        containing partner is the one cross-class operand, the principled
        exception the base algebra's ``(T, T) -> T`` cannot carry (#288).
        """
        partner: object = other  # narrowed alone, so ``other`` keeps the base algebra's type
        containing = type(partner).CONTAINS if isinstance(partner, BulkField) else ()
        if isinstance(self, containing):
            return NotImplemented
        return super().__add__(other)
'''

PLANTED_BULK_FIELD_CONTAINS = "    CONTAINS = ()\n"

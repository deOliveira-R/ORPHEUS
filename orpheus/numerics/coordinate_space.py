r"""Coordinate-space membership: the shape leg of the one parse.

An element of a space :math:`V` is handed to a verb as its coordinates, an
array. Before the verb computes anything, the array is parsed as an element
of the coordinate space :math:`\mathbb R^{\text{shape}}`. That parse has two
legs, which are two different facts:

* the **shape leg**, "``x`` is an element of :math:`\mathbb R^{\text{shape}}`",
  which needs no space, no axis and no metric, and lives here;
* the **identity leg**, "which space", the content comparison of a typed
  member with a :class:`~orpheus.numerics.space.FunctionSpace`, which lives
  on the space (:meth:`~orpheus.numerics.space.FunctionSpace.coordinates_of`).

Why a parse and not numpy's arithmetic
======================================

NumPy's ``*`` is the Hadamard product on the broadcast closure of shapes, a
larger domain than :math:`V \times V`. Handed two shapes that differ, it
silently extends one operand along the missing axes, and computes one of:

* an outer product, both operands extended (ERR-122's 0.225 for 0.095);
* the implicit pullback :math:`x \otimes 1`, a different member from the one
  passed, and with the wrong unit on a modal axis, where the pullback is
  :math:`u \otimes x`;
* an array that is not a member of the codomain at all.

The two extensions the mathematics does admit are named elsewhere: the tensor
action :math:`A \otimes I` along the axes an operator does not act on, placed
by axis label when the operator is built, and the pullback
:math:`\pi^* = u \otimes \mathrm{id}`
(:class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`). Nothing else
extends, and this module is what refuses the rest, on every run (``-O``
included): the check is the parse, not a guard.

The two verbs
=============

* :func:`coordinates`, the exact parse, for a scalar verb (a pairing, a norm,
  a mint, an integral) and for an engine minted with its full domain shape.
* :func:`leading_coordinates`, the parse of an action :math:`A \otimes I_P`:
  the leading axes must be the domain's, and an unlabelled trailing payload
  :math:`P` rides through.

Both return the operand itself, never a copy, so no product or reduction
downstream changes. This module sits below every numerics concept (the first
row of the concept order) so that a measure, a basis and a metric block, which
own their domain shapes but cannot name a space, parse with it too.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = ["NotInCoordinateSpace", "coordinates", "leading_coordinates"]


class NotInCoordinateSpace(ValueError):
    r"""An array is not an element of the coordinate space a verb acts on.

    A wrong shape is a wrong value of the right type, so this is a
    :class:`ValueError`; a non-array operand, where one is refused, is a
    :class:`TypeError` raised by its owner.
    """


def coordinates(x: Any, shape: tuple[int, ...], *, owner: str) -> Any:
    r"""Return ``x`` itself if it is an element of :math:`\mathbb R^{\text{shape}}`, or refuse it.

    Parameters
    ----------
    x : array_like
        The operand's coordinates.
    shape : tuple of int
        The exact shape of the coordinate space.
    owner : str
        The verb that parses, named in the refusal (the raise sits one frame
        below it).

    Raises
    ------
    NotInCoordinateSpace
        If ``np.shape(x) != shape``.
    """
    if np.shape(x) != shape:
        raise NotInCoordinateSpace(
            f"{owner}: values.shape {np.shape(x)!r} does not match space.shape "
            f"{shape!r} (an element takes the exact shape)."
        )
    return x


def leading_coordinates(x: Any, shape: tuple[int, ...], *, owner: str) -> Any:
    r"""Return ``x`` itself if its leading axes are :math:`\mathbb R^{\text{shape}}`, or refuse it.

    The parse of an action :math:`A \otimes I_P`: ``A`` acts on the leading
    axes and an unlabelled trailing payload :math:`P` rides through.

    Raises
    ------
    NotInCoordinateSpace
        If ``np.shape(x)[:len(shape)] != shape``.
    """
    if np.shape(x)[: len(shape)] != shape:
        raise NotInCoordinateSpace(
            f"{owner}: an array of shape {np.shape(x)!r} does not have the leading "
            f"shape {shape!r} that an action acts on (A ⊗ I_P)."
        )
    return x

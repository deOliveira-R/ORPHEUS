"""The boundary-condition tag: a solver-agnostic declaration attached to a geometry's surfaces.

A :class:`BC` is ``(kind, params)``: a name each method resolves through its
own registry, plus the numeric parameters a law needs. It belongs to the
geometry layer because boundaries are defined at the geometry; the typed laws
it resolves to (:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` and its
realizations) live beside it in this package.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar

__all__ = ["BC"]


@dataclass(frozen=True)
class BC:
    """Solver-agnostic boundary condition declaration.

    A lightweight tag attached to geometry surfaces.  The geometry module
    makes no assumptions about what a given ``kind`` means — semantics
    are resolved by each solver's augmented mesh at construction time
    via its ``BC_REGISTRY``.

    Parameters
    ----------
    kind : str
        Boundary condition identifier (e.g. ``"vacuum"``,
        ``"reflective"``, ``"white"``).  Each solver defines which
        kinds it supports.
    params : dict[str, float]
        Optional numeric parameters (e.g. ``{"albedo": 0.7}``).
    """

    kind: str
    params: dict[str, float] = field(default_factory=dict)

    # ── the three parameter-free tags, as named constants ────────────────
    #
    # Declared here and BOUND after the class body (a class cannot
    # instantiate itself inside its own body).  ``ClassVar`` is what keeps
    # them out of the dataclass field list — and what makes them real
    # documented attributes: an annotation in the class body lands in
    # ``__annotations__``, so autodoc emits a ``py:attribute`` for each and
    # ``:obj:`BC.vacuum <orpheus.geometry.mesh.BC.vacuum>``` resolves.
    #
    # ⛔ Until 2026-08-10 these existed ONLY as post-class assignments
    # carrying ``# type: ignore[attr-defined]``.  Every consequence of that
    # was invisible until an instrument looked: the type checker had to be
    # silenced three times, autodoc had nothing to document, and every
    # cross-reference to them — 7 sites over 3 pages — rendered as plain
    # text.  #346 W1 (`52650a86`) then QUALIFIED those references to the
    # full dotted path, which is the right spelling and made the graph
    # report them DEAD: a bare ``BC.vacuum`` had been resolvable by
    # Sphinx's suffix search, a qualified one needs a real object.  The
    # honest fix is not to un-qualify the references; it is to make the
    # attribute exist to the tooling, which is what these three lines do.
    vacuum: ClassVar["BC"]
    """The no-return tag :math:`\\alpha = 0` — nothing enters from outside."""
    reflective: ClassVar["BC"]
    """The specular tag :math:`\\alpha = 1` — every ordinate mirrors back."""
    white: ClassVar["BC"]
    """The isotropic-return tag — the re-emission closure, not a mirror."""

    def __repr__(self) -> str:
        if self.params:
            return f"BC({self.kind!r}, {self.params!r})"
        return f"BC({self.kind!r})"

    def to_alpha(self) -> float:
        r"""Map this BC tag to a continuous specular albedo
        :math:`\alpha \in [0, 1]`.

        The trajectory_resolvent / Birkhoff–Sinai billiard family
        parametrises specular boundary conditions on a continuous
        albedo: :math:`\alpha = 0` is vacuum (no return), :math:`\alpha
        = 1` is perfect specular reflection, :math:`\alpha \in (0, 1)`
        is partial reflection. This method translates the production
        :class:`BC` tag-system to that scalar:

        * :data:`BC.vacuum` → ``0.0``
        * :data:`BC.reflective` → ``1.0``
        * ``BC("partial", {"albedo": x})`` → ``x``

        Other tags (``"white"``, Marshak diffuse, …) raise
        :class:`NotImplementedError` — they require a different
        closure structure that the specular-albedo parametrisation
        does not represent.

        Returns
        -------
        float
            The continuous-albedo equivalent of this BC tag.

        Raises
        ------
        ValueError
            If the BC kind is ``"partial"`` and the ``"albedo"`` key
            is missing from :attr:`params`.
        NotImplementedError
            If the BC kind has no specular-albedo equivalent
            (e.g. ``"white"``).
        """
        if self.kind == "vacuum":
            return 0.0
        if self.kind == "reflective":
            return 1.0
        if self.kind == "partial":
            try:
                return float(self.params["albedo"])
            except KeyError as exc:
                raise ValueError(
                    f"BC.to_alpha: BC('partial', ...) is missing the "
                    f"'albedo' parameter; got params={self.params!r}"
                ) from exc
        raise NotImplementedError(
            f"BC.to_alpha: BC kind {self.kind!r} has no specular-albedo "
            f"equivalent. Tags supported today: vacuum, reflective, partial."
        )


# The bindings for the three ClassVars declared in the class body — see the
# note there for why the declaration and the binding are separated, and for
# what the three retired ``# type: ignore[attr-defined]`` comments were
# hiding.
BC.vacuum = BC("vacuum")
BC.reflective = BC("reflective")
BC.white = BC("white")

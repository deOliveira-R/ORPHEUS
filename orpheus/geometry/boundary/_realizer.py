r"""BoundaryRealizer Protocol + BoundaryRealizerRegistry.

Per Grand Report v3 §16A.3 lines 2841–2860, the third layer of the
boundary architecture (after trace structure + physical law) is the
**method realisation**: a per-transport-method strategy that turns a
BC descriptor (:class:`BoundaryTraceLaw`) into a
:class:`~orpheus.numerics.operator.LinearOperator` consumable by that
method's sweep / solver.

The three-layer split matters because the same physical law (vacuum,
specular reflection, white, …) is realized by **structurally
different** linear operators in each transport method:

* SN realizes vacuum as a sparse per-ordinate
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor` on
  the inflow ordinates of the affected face;
* MoC realizes the same vacuum law as zeroing the entering track
  boundary fluxes;
* MC realizes it by killing particles that exit;
* CP realizes it as zero rows in the boundary-to-region coupling
  matrix; …

The :class:`BoundaryRealizer` Protocol lives in this module; the
concrete realizers ship in per-method subpackages:

* :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` —
  functional realizer for SN (dispatches by ``isinstance(law, ...)``
  to the Wave-0 / Wave-1 primitives).
* :class:`~orpheus.moc.boundary_realizer.MoCBoundaryRealizer`,
  :class:`~orpheus.mc.boundary_realizer.MCBoundaryRealizer`,
  :class:`~orpheus.cp.boundary_realizer.CPBoundaryRealizer`,
  :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
  — stub realizers (``NotImplementedError`` with grep-able marker)
  registered for the day each method adopts the unified BC
  architecture.

The :class:`BoundaryRealizerRegistry` is a stand-alone registry (not
mounted on :class:`~orpheus.numerics.registry.RegistryMixin`) because
realizers are NOT a hierarchy of subclasses — they're independent
strategies keyed by ``method_name``. The registry's collision rule
fails LOUDLY at import time (``BoundaryRealizerRegistryError``), so
duplicate-registration bugs surface during package import rather
than at first ``realize()`` call.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 5 — C5.1 / C5.2.
* Grand Report v3 §16A.4 lines 2864–2876 (``realize(law, method_space)``
  vocabulary).
* Grand Report v3 §16A.11 lines 3252–3257 (the two-registry design:
  law-registry on :class:`BoundaryTraceLaw`, realizer-registry here).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator


__all__ = [
    "BoundaryRealizer",
    "BoundaryRealizerRegistry",
    "BoundaryRealizerRegistryError",
]


class BoundaryRealizerRegistryError(KeyError):
    """Raised on :class:`BoundaryRealizerRegistry` collisions / lookup misses.

    Subclasses :class:`KeyError` (the natural Python type for
    registry-key issues). The collision case fires at import time
    when two realizers are registered under the same
    ``method_name``; the miss case fires when
    :meth:`BoundaryRealizerRegistry.get` is called with a key that
    was never registered (typically meaning the consumer didn't
    import the realizer module yet).
    """


@runtime_checkable
class BoundaryRealizer(Protocol):
    r"""Method-specific realisation of a boundary law.

    Implementors live in per-method subpackages
    (``orpheus.sn.boundary_realizer``, etc.) and self-register via
    :meth:`BoundaryRealizerRegistry.register` at import time.

    The :meth:`realize` method takes a method-agnostic boundary law
    (the BC descriptor :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`)
    and a method-specific *method space* — a lightweight container
    holding whatever discretisation metadata the realizer needs
    (quadrature, mesh, trace masks, …) — and returns a
    :class:`~orpheus.numerics.operator.LinearOperator` whose 1-arg
    :meth:`apply(psi)` matches the legacy 2-arg
    :meth:`law.apply(psi, quadrature)` semantics.

    The intent of the Protocol is **structural** (the dispatch
    rules + the registry registration), not algorithmic — every
    realizer's :meth:`realize` body is method-specific. Wave 5 ships
    the SN realizer as the functional reference; the other methods
    ship stubs (NotImplementedError) that hold the architecture
    end-to-end so Wave 7 can rebase the legacy concretes onto
    :class:`BoundaryTraceLaw` without leaving holes in the dispatch
    table.

    Attributes
    ----------
    method_name : str
        Stable identifier (``"SN"``, ``"MoC"``, ``"MC"``, ``"CP"``,
        ``"diffusion"``, …) used as the registry key.
    """

    method_name: str

    def realize(
        self,
        law: Any,
        method_space: Any,
    ) -> "LinearOperator":
        """Return a 1-arg :class:`LinearOperator` realising ``law`` for this method."""
        ...


class BoundaryRealizerRegistry:
    """Method-keyed registry of :class:`BoundaryRealizer` strategies.

    Used as ``BoundaryRealizerRegistry.register("SN")(MyRealizer)``
    to add an entry, ``BoundaryRealizerRegistry.get("SN")`` to look
    one up, and ``BoundaryRealizerRegistry.method_names()`` to
    enumerate every registered method. Collisions raise
    :class:`BoundaryRealizerRegistryError` at import time —
    surfacing the double-registration loudly rather than letting
    one realizer silently shadow another.

    The registry is class-level state (singleton per process); a
    process-local design suffices because realizers self-register at
    module import time and the registry is never cleared during
    normal operation. A future test-harness reset hook (if needed)
    can clear the ``_registry`` dict directly.
    """

    _registry: dict[str, type[BoundaryRealizer]] = {}

    @classmethod
    def register(cls, method_name: str):
        """Class decorator: register the realizer under ``method_name``.

        Collisions raise :class:`BoundaryRealizerRegistryError`
        immediately — duplicate registrations are bugs, not
        overrides. If a consumer needs to override a registration,
        they should clear the entry explicitly via ``del
        BoundaryRealizerRegistry._registry[method_name]`` before
        registering the replacement.
        """

        def decorator(
            realizer_cls: type[BoundaryRealizer],
        ) -> type[BoundaryRealizer]:
            if method_name in cls._registry:
                raise BoundaryRealizerRegistryError(
                    f"BoundaryRealizerRegistry already contains an entry "
                    f"for method_name={method_name!r}: "
                    f"{cls._registry[method_name].__name__}. Refusing to "
                    f"overwrite with {realizer_cls.__name__}."
                )
            cls._registry[method_name] = realizer_cls
            return realizer_cls

        return decorator

    @classmethod
    def get(cls, method_name: str) -> type[BoundaryRealizer]:
        """Return the realizer class for ``method_name``.

        Raises :class:`BoundaryRealizerRegistryError` if no
        realizer is registered under that key — typically meaning
        the consumer didn't import the realizer module yet (the
        stubs self-register on import, so the fix is usually a
        missing ``import orpheus.<method>``).
        """
        try:
            return cls._registry[method_name]
        except KeyError as exc:
            available = sorted(cls._registry.keys())
            raise BoundaryRealizerRegistryError(
                f"No BoundaryRealizer registered for "
                f"method_name={method_name!r}. Available: {available}."
            ) from exc

    @classmethod
    def method_names(cls) -> tuple[str, ...]:
        """All registered method names, sorted lexicographically."""
        return tuple(sorted(cls._registry.keys()))

r"""The labels a test gives a hand-built coupled field: ``s0``, ``s1``, ... one per member system. A helper RETURNS data and asserts nothing."""

from __future__ import annotations

from collections.abc import Sequence


def labels_for(systems: "Sequence[object] | int") -> tuple[str, ...]:
    """One summand label ``s<i>`` per member: ``systems`` is the members' sequence or their count (a count spares a call site from evaluating its members twice)."""
    count = systems if isinstance(systems, int) else len(systems)
    return tuple(f"s{i}" for i in range(count))


def hub_labels(systems: "Sequence[object] | int") -> tuple[str, ...]:
    """The summand labels of the production coupled sum for ``systems`` members: the builder's own constants, System A first, System B second."""
    from orpheus.sn.coupled_system import SYSTEM_A_LABEL, SYSTEM_B_LABEL

    count = systems if isinstance(systems, int) else len(systems)
    if not 1 <= count <= 2:
        raise ValueError(f"the production coupled sum has one or two summands; asked for {count}")
    return (SYSTEM_A_LABEL, SYSTEM_B_LABEL)[:count]

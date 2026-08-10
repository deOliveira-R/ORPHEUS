r"""The registry NAME of a Peierls continuous reference — one rule, one home.

Why this module exists (#345)
=============================

A Peierls reference's registry name is determined entirely by its
*identity*: which shape, how many energy groups, how many spatial regions,
and — for a hollow curvilinear case — the **dimensionless cavity ratio**
:math:`r_0/R`. Nothing else. So the name is a *function*, and before #345 it
was written out as **14 separate expressions across 4 modules** (`[M]`
2026-08-09): 5 authoritative ``ContinuousReferenceSolution`` stamps
(``cylinder`` ×2, ``sphere`` ×2, ``slab``), 1 more in ``cases``' unified slab
path, 3 lazy-builder keys (one a bare literal), and 5 capability-matrix row
templates (one likewise a literal). The :math:`r_0` tag inside them had **6**
spellings of its own.

`[M]` — **two of those tag spellings already disagreed**, and only a
coincidence hid it:

======================================  =========================  ==========
site                                    tag formula                correct?
======================================  =========================  ==========
``cylinder.py`` / ``sphere.py`` stamp   ``round(r0_over_R * 100)``  ✅ (×2)
``cases.continuous_case_builders``      ``round(r0/R_out * 100)``   ✅ (×2)
``cases.capability_rows``               ``round(r0 * 100)``         ⛔ (×2) no ``/R_out``
======================================  =========================  ==========

They agree today because every outer radius in
:data:`~orpheus.derivations.continuous.flat_source_cp.cylinder._RADII` is
``1.0`` — for *every* region count, which is exactly what makes the defect
survivable: a contributor adding a non-unit outer radius would silently
publish a capability matrix naming references that are not registered, with
no gate anywhere able to see it (the join is by name, and the name would be
wrong on both sides of nothing).

The deeper defect the tag divergence was a symptom of
-----------------------------------------------------

The sweep is over a **dimensionless ratio** — every docstring in the family
says *"hollow cylinder F.4 at* :math:`r_0/R \in \{0.1, 0.2, 0.3\}`\ *"* — but
it was carried in ``build_two_surface_case(inner_radius=...)``, a parameter
the signature declares to be an absolute **length**, which then divided by
``R_out`` to recover the ratio. Passing ``0.1`` therefore meant
:math:`r_0 = 0.1\,\text{cm}` to the API and :math:`r_0/R = 0.1` to the
caller, and the two coincided only at :math:`R = 1`. That is
`coding-elegance` Pattern 3 in its sharpest form — **an intermediate whose
units are wrong is not a naming problem, it is a physics problem** — and
Pattern 7: a convention re-applied at every consumer instead of fixed once
at the definition site.

#345 fixed the units at the definition site (``build_two_surface_case`` now
takes ``r0_over_R``, the quantity the domain actually names, and performs no
division) and hoisted the enumeration to a single grid. This module holds
the naming half.

What this buys
--------------

The drift the capability matrix's promised cross-check was meant to *detect*
is now **unspellable**: there is one enumeration and one name rule, so a row
and a registry key cannot disagree. Per `coding-standards`' rewire-demotion
clause, the surviving gate in
``tests/derivations/test_capability_matrices.py`` says so plainly rather than
keeping an authoritative name for a comparison that can no longer fail.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["SHAPE_TAGS", "ShippedReference", "reference_name"]

#: Registry-name fragment per shape. The keys are the ``shape`` vocabulary
#: :func:`~orpheus.derivations.continuous.peierls_nystrom.cases.build_two_surface_case`
#: accepts; the values are the abbreviations the registry names carry.
SHAPE_TAGS: dict[str, str] = {
    "slab": "slab",
    "cylinder-1d": "cyl1D",
    "sphere-1d": "sph1D",
}


def reference_name(
    shape: str,
    *,
    n_groups: int,
    n_regions: int,
    r0_over_R: float | None = None,
) -> str:
    r"""The canonical registry name of a Peierls continuous reference.

    Parameters
    ----------
    shape
        One of :data:`SHAPE_TAGS` — ``"slab"``, ``"cylinder-1d"`` or
        ``"sphere-1d"``.
    n_groups
        Energy-group count.
    n_regions
        Spatial-region count.
    r0_over_R
        The **dimensionless** cavity ratio :math:`r_0/R` for a hollow
        curvilinear case, or ``None`` for a slab / solid case. It is a
        ratio, never a length: passing an absolute :math:`r_0` produces a
        name that is right only when :math:`R = 1`, which is the exact
        defect #345 removed.

    Returns
    -------
    str
        ``peierls_<shape>[_hollow]_<ng>eg_<nr>rg[_r0_<tag>]`` with
        ``tag = round(100 * r0_over_R)`` zero-padded to two digits.

    Raises
    ------
    ValueError
        On an unknown ``shape``, or an ``r0_over_R`` outside ``(0, 1)`` —
        a cavity that is not strictly inside its own outer surface is not a
        hollow case, and the most likely way to produce one is to pass a
        length where the ratio belongs.
    """
    try:
        tag = SHAPE_TAGS[shape]
    except KeyError:
        raise ValueError(
            f"reference_name: unknown shape {shape!r}; "
            f"expected one of {sorted(SHAPE_TAGS)}"
        ) from None

    stem = f"peierls_{tag}"
    if r0_over_R is None:
        return f"{stem}_{n_groups}eg_{n_regions}rg"

    if not 0.0 < r0_over_R < 1.0:
        raise ValueError(
            f"reference_name: r0_over_R must be a RATIO in (0, 1), got "
            f"{r0_over_R!r}. A hollow case's cavity lies strictly inside "
            f"its outer surface; a value at or above 1 is almost always an "
            f"absolute radius passed where the dimensionless r_0/R belongs."
        )
    return f"{stem}_hollow_{n_groups}eg_{n_regions}rg_r0_{round(r0_over_R * 100):02d}"


@dataclass(frozen=True, slots=True)
class ShippedReference:
    """One shipped Peierls reference, by IDENTITY alone.

    These four fields are exactly what determines both the registry name and
    the construction call. Everything else a consumer needs — the accuracy
    prose, the closure label, the production-status string of the capability
    matrix — legitimately differs per consumer and stays at the consumer
    (`[[lessons-L32]]`: factor out the invariant slice, leave the varying one
    where it varies).
    """

    shape: str
    n_groups: int
    n_regions: int
    r0_over_R: float | None = None

    @property
    def name(self) -> str:
        """The registry name — the single spelling, via :func:`reference_name`."""
        return reference_name(
            self.shape,
            n_groups=self.n_groups,
            n_regions=self.n_regions,
            r0_over_R=self.r0_over_R,
        )

    @property
    def ng_key(self) -> str:
        """The ``"1g"`` / ``"2g"`` cross-section-library key."""
        return f"{self.n_groups}g"

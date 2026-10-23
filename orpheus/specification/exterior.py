r"""The reference's wall value ``WallAmplitudes(specular, diffuse)``.

The reference specification states a wall by the coordinates :math:`(s, d)` of
its linear part in the triangle K of sub-Markov walls (the positive cone cut by
:math:`s + d \le 1`; its vertices are the vacuum, the mirror and the white wall): the fraction :math:`s` of the incoming
current returned specularly and the fraction :math:`d` returned diffusely, as
plain data (the interface between references and production is plain data).
A wall that returns more than it receives is not physical, so both are in
:math:`[0, 1]` with :math:`s + d \le 1`; a value outside that is a parse
refusal (:class:`ValueError`, :class:`TypeError` for a non-real).

The vocabulary is closed: a wall is vacuum :math:`(0, 0)`, specular
:math:`(s, 0)` or diffuse :math:`(0, d)`. A combined wall (both coordinates
positive) is a scope boundary that issue #650 lifts by widening the
vocabulary, so it raises :class:`NotImplementedError`, not :class:`ValueError`:
no better input lifts it. The bounds are checked first, so an unphysical
combined pair such as :math:`(0.6, 0.6)` is a :class:`ValueError`.
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.values.content import ContentIdentity
from orpheus.values.scalars import parse_real


@dataclass(frozen=True, eq=False)
class WallAmplitudes(ContentIdentity):
    """The linear part of a wall as coordinates in the triangle K: ``specular`` and ``diffuse`` returned fractions."""

    specular: float
    diffuse: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "specular", parse_real(self.specular, "WallAmplitudes: the specular amplitude"))
        object.__setattr__(self, "diffuse", parse_real(self.diffuse, "WallAmplitudes: the diffuse amplitude"))
        # The simplex: s >= 0, d >= 0, s + d <= 1 (the per-coordinate upper bounds follow from these).
        if not (self.specular >= 0.0 and self.diffuse >= 0.0 and self.specular + self.diffuse <= 1.0):
            raise ValueError(
                f"WallAmplitudes(specular={self.specular}, diffuse={self.diffuse}): a returned amplitude outside "
                "[0, 1], or returning more than it receives, is not a physical wall"
            )
        if self.specular > 0.0 and self.diffuse > 0.0:
            # SCOPE-BOUNDARY[guard] -- machinery: the reference vocabulary widened to the two-chain wall
            # s*Specular + d*Diffuse (plan 2.9, from_amplitudes). ruling: D3, DRAFT 3.1, approved by the
            # user 2026-10-23. revisit: #650 (state prescribed inflow, combined walls and zero flux).
            raise NotImplementedError(
                f"WallAmplitudes(specular={self.specular}, diffuse={self.diffuse}): a wall is specular or diffuse "
                "in the closed vocabulary; the combined wall waits on #650"
            )

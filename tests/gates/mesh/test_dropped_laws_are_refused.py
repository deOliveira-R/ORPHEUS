"""A lone translation face is not a slab: the geometry's pairing rule, the half of MC's old refusal rows that is no method's.

Until the wave-2 carve MC refused "a periodic right face with a non-periodic left law" on a slab. A translation deck needs its
partner (an involution), so the geometry refuses ``TranslationDeck`` against a tag with ``DeckError`` and MC never sees the dropped
left law the guard was written for; the row below pins that.

This file held the dropped-law refusals of CP, MoC and MC (#513, #514; the P1 specification's S3.12 to S3.14) until the layer-DAG
re-filing (Phase 1b, 2026-10-26). Each refusal is one method's claim, so each class moved to its method's gates, bodies unchanged:
``tests/gates/cp/test_cp_dropped_laws_are_refused.py``, ``tests/gates/moc/test_moc_dropped_laws_are_refused.py`` and
``tests/gates/mc/test_mc_dropped_laws_are_refused.py``, with their shared mesh builders in
``tests/gates/transport/_dropped_laws_helpers.py``. The row here executes only the geometry, and the class keeps its name so its node
id is unchanged.
"""
from __future__ import annotations

import pytest

from orpheus.geometry import StructuredGeometry, TranslationDeck
from orpheus.geometry.deck import DeckError

pytestmark = pytest.mark.foundation


class TestMC:
    def test_a_lone_translation_face_is_not_a_slab(self):
        """The old row "a periodic right face with a non-periodic left law" cannot be declared: a translation deck needs its partner (the
        geometry refuses ``TranslationDeck`` against a tag, ``DeckError``), so MC never sees the dropped left law the guard was written for."""
        with pytest.raises(DeckError):
            StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=TranslationDeck())

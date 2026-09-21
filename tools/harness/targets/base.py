"""What the pipeline needs from a harness.

A harness is a target runtime: it decides, per source kind, where a page
lands, what shape it takes there, and whether it is loaded into every
session. One implementation exists today; the Protocol is the seam the user
required for the next one (ruling 2026-09-20), and its rent is
the enforced boundary — the neutral modules import nothing below this file
and name no harness, which ``tests/tools/test_harness_generator.py`` asserts.
"""
from __future__ import annotations

from pathlib import Path
from typing import Protocol

from ..render import Rendered
from ..source import Kind, Page


class Harness(Protocol):
    @property
    def name(self) -> str: ...
    @property
    def kinds(self) -> frozenset[Kind]: ...      # the kinds this harness realises; the pipeline skips the others
    @property
    def roots(self) -> tuple[Path, ...]: ...     # the directories and files it owns: where orphans are looked for

    def target(self, page: Page) -> Path: ...
    def render(self, page: Page, body: str) -> Rendered: ...   # body: the page's body, links already re-pointed at the target
    def always_on(self, page: Page) -> bool: ...              # loaded into every session and every inheriting dispatch
    def installed_always_on(self) -> tuple[Path, ...]: ...    # always-on files another tool installed under the roots; no page produces them

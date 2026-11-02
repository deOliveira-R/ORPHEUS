"""The memo's pins as the store keeps them: one row form (the memo's own) and paths relative to the repository.

A pin is the memo's own object (:data:`orpheus.numerics.traced_memo.Pin`); its row is the memo's
(``Manifest.to_json``: ``[kind, *fields]``) and its decoder the memo's (``Manifest.from_json``). The one
difference between a pin as the memo makes it and as the store keeps it is where a path points: the memo's
``DefPin`` and ``ModulePin`` name a file relative to a ``sys.path`` entry (or absolutely), and its ``DataPin``,
``PresencePin`` and ``ListingPin`` absolutely; the store keeps all of them relative to the repository, so a store
survives a move of the checkout. :func:`relocated` is the one function that changes a path.
"""
from __future__ import annotations

import json
from collections.abc import Callable

from orpheus.numerics import traced_memo as memo

#: Pins that name a file's bytes: judged by the file's identity first (stage one), by the pin second.
CONTENT_PINS = (memo.DefPin, memo.ModulePin, memo.DataPin)
#: Every pin that names a path.
PATH_PINS = (*CONTENT_PINS, memo.PresencePin, memo.ListingPin)
PathPin = memo.DefPin | memo.ModulePin | memo.DataPin | memo.PresencePin | memo.ListingPin


def file_of(pin: memo.Pin) -> str | None:
    """The file a content pin judges, or ``None`` for any other pin."""
    return pin.path if isinstance(pin, CONTENT_PINS) else None


def relocated(pin: memo.Pin, path: Callable[[PathPin], str | None]) -> memo.Pin | None:
    """``pin`` with its path replaced by ``path(pin)``; ``None`` when ``path`` returns ``None`` (the pin is
    dropped). A pin that names no path passes through."""
    if not isinstance(pin, PATH_PINS):
        return pin
    moved = path(pin)
    return None if moved is None else pin._replace(path=moved)


def encode(pin: memo.Pin) -> str:
    return json.dumps(memo.Manifest(frozenset({pin})).to_json()[0])


def decode(row: str) -> memo.Pin:
    (pin,) = memo.Manifest.from_json([json.loads(row)]).pins
    return pin

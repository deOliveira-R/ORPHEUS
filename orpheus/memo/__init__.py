"""The traced memo: a pure function of content, memoised on disk, with the import tracer that pins what it ran.

Layer law: this package imports only the standard library, numpy and :mod:`orpheus.values`; nothing below it (``geometry``,
``values``) imports it. It holds no import of its own, so that ``import orpheus.memo.traced_memo`` loads the memo and its
package only.
"""

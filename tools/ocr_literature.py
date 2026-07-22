#!/usr/bin/env python3
"""Mistral OCR sidecars for the literature library (`scratch/literature/`).

Purpose
-------
Give the literature-researcher a *text* surface for the scanned library so
extraction becomes selection-from-input instead of generation-as-output (the
2026-07-22 content-filter incident: transcribing a textless scan as model
output trips the API's output classifier; text read from disk does not).

The SCAN is the SSOT. Sidecars are for search, navigation, and prose
quotation ONLY — every load-bearing equation is spot-verified against the
rendered page (Read with ``pages=``) before it enters a theory page or a
solver (the ERR-032-class hazard: plausible-but-wrong reference equations).

Provenance
----------
OCR plumbing lifted from ``~/git/RX7_workshop_manual/manual_figures/extract.py``
(client init :62, model :76, ``ocr.process`` :85, cache-idempotent
``ocr_cached`` :91). The OpenCV figure-cropping stages are not needed here.

Run (Host env — always ``.venv/bin/python``)
--------------------------------------------
  .venv/bin/python tools/ocr_literature.py --validate --glob 'Sanchez(1982)*'
  .venv/bin/python tools/ocr_literature.py --glob 'Sanchez(1982)*'
  .venv/bin/python tools/ocr_literature.py                # full library sweep

Behavior
--------
- Raw responses cache to ``scratch/literature_ocr/<stem>.mocr.json`` —
  re-runs never re-bill (~$1 per 1000 pages). Sidecars
  (``<stem>.md``) are regenerated from cache on every run, so emitter
  improvements are free.
- ``## p. <N>`` headings are 1-based PDF pages — the same numbers the Read
  tool's ``pages=`` parameter takes. Page boundaries are never flattened
  (corpus citation discipline is page-anchored).
- Printed-page mapping (e.g. journal p. 481 = PDF p. 2) is detected from
  OCR'd headers/footers: only the true printed label co-varies linearly
  with the PDF index, so the modal offset filters volume/year noise.
- ``MISTRAL_API_KEY`` loads from the ORPHEUS-root ``.env`` (gitignored).
"""

from __future__ import annotations

import argparse
import base64
import fnmatch
import json
import os
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "scratch" / "literature"
DEFAULT_OUT = REPO_ROOT / "scratch" / "literature_ocr"
COST_PER_PAGE_USD = 1.0 / 1000.0  # mistral-ocr-latest, 2026-07 pricing

# ----------------------------- Mistral client -----------------------------


def load_env_key() -> None:
    """Populate MISTRAL_API_KEY from the ORPHEUS-root .env (stdlib parse,
    RX7 pattern). Existing environment wins; the file is never echoed."""
    env_path = REPO_ROOT / ".env"
    if not env_path.is_file():
        return
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        os.environ.setdefault(key.strip(), val.strip().strip('"').strip("'"))


def get_client():
    load_env_key()
    if not os.environ.get("MISTRAL_API_KEY"):
        sys.exit("MISTRAL_API_KEY not set (checked environment and the repo-root .env)")
    from mistralai.client import Mistral  # SDK 2.x home (1.x exported top-level)

    return Mistral(api_key=os.environ["MISTRAL_API_KEY"])


def run_ocr(client, pdf_path: Path, pages: list[int] | None = None) -> dict:
    """Document OCR via base64 data-URL. Headers/footers are extracted as
    separate per-page fields — that is where printed page numbers live."""
    b64 = base64.b64encode(pdf_path.read_bytes()).decode("utf-8")
    kwargs: dict[str, object] = dict(
        document={
            "type": "document_url",
            "document_url": f"data:application/pdf;base64,{b64}",
        },
        model="mistral-ocr-latest",
        include_image_base64=False,
        include_blocks=False,
        table_format="markdown",
        extract_header=True,
        extract_footer=True,
    )
    if pages is not None:
        kwargs["pages"] = pages
    resp = client.ocr.process(**kwargs)
    if hasattr(resp, "model_dump"):
        return resp.model_dump()
    return json.loads(resp.json())


def ocr_cached(
    client, pdf_path: Path, cache_path: Path, force: bool = False
) -> tuple[dict, bool]:
    """Full-document OCR through the cache. Returns (data, billed)."""
    if not force and cache_path.exists():
        return json.loads(cache_path.read_text()), False
    data = run_ocr(client, pdf_path)
    cache_path.write_text(json.dumps(data))
    return data, True


# ----------------------------- telemetry -----------------------------


def text_layer_chars(pdf_path: Path, n_pages: int = 3) -> int:
    """Extracted-text length over the first pages (pypdf). Telemetry only —
    OCR runs on everything per the 2026-07-22 ruling. Returns -1 on failure."""
    try:
        from pypdf import PdfReader

        reader = PdfReader(str(pdf_path))
        return sum(
            len((page.extract_text() or "").strip())
            for page in reader.pages[:n_pages]
        )
    except Exception:
        return -1


_YEAR_RANGE = range(1900, 2100)


def printed_page_offset(pages: list[dict]) -> tuple[int, int] | None:
    """Detect ``printed = pdf_1based + offset`` from headers/footers.

    Volume numbers, years, and page-range citations appear on every page at
    a CONSTANT value, so their offset scatters as the PDF index advances;
    only the true printed label tracks the index linearly and produces a
    constant offset. The modal offset therefore filters the noise. Years
    (1900–2099) are dropped outright. Returns (offset, n_supporting) when
    the mode has >=5 supporters covering >=50% of number-bearing pages.
    """
    offsets: Counter[int] = Counter()
    pages_with_numbers = 0
    for page in pages:
        pdf_no = page.get("index", 0) + 1
        text = f"{page.get('header') or ''} {page.get('footer') or ''}"
        candidates = {
            int(tok)
            for tok in re.findall(r"\b(\d{1,4})\b", text)
            if int(tok) not in _YEAR_RANGE
        }
        if candidates:
            pages_with_numbers += 1
            for printed in candidates:
                offsets[printed - pdf_no] += 1
    if not offsets:
        return None
    offset, support = offsets.most_common(1)[0]
    if support >= 5 and support >= 0.5 * pages_with_numbers:
        return offset, support
    return None


# ----------------------------- sidecar emitter -----------------------------


def write_sidecar(
    out_path: Path, pdf_path: Path, data: dict, chars: int, ocr_day: str
) -> int:
    pages = sorted(data.get("pages") or [], key=lambda p: p.get("index", 0))
    mapping = printed_page_offset(pages)
    if mapping is not None:
        offset, support = mapping
        example_pdf = pages[len(pages) // 2].get("index", 0) + 1
        map_line = (
            f"printed p. ≈ PDF p. + {offset} "
            f"(heuristic, {support} supporting pages; e.g. PDF p. {example_pdf} "
            f"= printed p. {example_pdf + offset}) — verify before citing"
        )
    else:
        map_line = (
            "none detected — cite by PDF page, or read the printed number "
            "off the page's header/footer lines below"
        )
    if chars < 0:
        layer = "unknown (pypdf failed)"
    elif chars < 100:
        layer = f"none — a scan ({chars} chars over the first 3 pages)"
    else:
        layer = f"yes ({chars} chars over the first 3 pages)"

    rel = pdf_path.relative_to(REPO_ROOT)
    lines = [
        f"# OCR sidecar — {pdf_path.stem}",
        "",
        f"> **Source (SSOT):** `{rel}` — the scan is the single source of",
        "> truth. This sidecar is for search, navigation, and prose quotation",
        "> ONLY. **Spot-verify every load-bearing equation against the rendered",
        "> page** (Read with `pages=`) before it enters a theory page or a",
        "> solver (ERR-032-class hazard).",
        ">",
        f"> OCR: Mistral `mistral-ocr-latest`, {ocr_day}. Text layer in source: {layer}.",
        "> Page headings `## p. N` are 1-based PDF pages — the numbers Read's",
        "> `pages=` takes.",
        f"> Printed-page mapping: {map_line}.",
        "",
    ]
    for page in pages:
        pdf_no = page.get("index", 0) + 1
        lines.append(f"## p. {pdf_no}")
        lines.append("")
        header = (page.get("header") or "").strip()
        footer = (page.get("footer") or "").strip()
        if header:
            lines.append(f"*[printed header: {header}]*")
            lines.append("")
        lines.append((page.get("markdown") or "").rstrip())
        lines.append("")
        if footer:
            lines.append(f"*[printed footer: {footer}]*")
            lines.append("")
    out_path.write_text("\n".join(lines))
    return len(pages)


# ----------------------------- driver -----------------------------


def select_pdfs(input_dir: Path, pattern: str) -> list[Path]:
    return sorted(
        p
        for p in input_dir.iterdir()
        if p.is_file() and fnmatch.fnmatchcase(p.name.lower(), pattern.lower())
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Mistral OCR sidecars for the literature library"
    )
    ap.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--glob", default="*.pdf", help="filename pattern (case-insensitive)")
    ap.add_argument(
        "--validate",
        action="store_true",
        help="OCR page 0 only, print it, write nothing (cheap pipeline check)",
    )
    ap.add_argument("--force", action="store_true", help="re-OCR even when cached")
    args = ap.parse_args()

    pdfs = select_pdfs(args.input_dir, args.glob)
    if not pdfs:
        sys.exit(f"no PDFs match {args.glob!r} in {args.input_dir}")
    client = get_client()

    if args.validate:
        for pdf in pdfs:
            print(f"=== validate (page 0 only): {pdf.name} ===")
            data = run_ocr(client, pdf, pages=[0])
            for page in data.get("pages") or []:
                md = page.get("markdown") or ""
                print(f"--- header: {page.get('header')!r}")
                print(f"--- footer: {page.get('footer')!r}")
                print(f"--- markdown ({len(md)} chars) ---")
                print(md)
        return 0

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ocr_day = date.today().isoformat()
    total_pages = billed_pages = 0
    failures: list[tuple[str, str]] = []
    for i, pdf in enumerate(pdfs, 1):
        chars = text_layer_chars(pdf)
        cache = args.out_dir / f"{pdf.stem}.mocr.json"
        sidecar = args.out_dir / f"{pdf.stem}.md"
        try:
            data, billed = ocr_cached(client, pdf, cache, force=args.force)
            n = write_sidecar(sidecar, pdf, data, chars, ocr_day)
        except Exception as exc:  # one bad PDF must not kill the sweep
            failures.append((pdf.name, f"{type(exc).__name__}: {exc}"))
            print(f"[{i}/{len(pdfs)}] {pdf.name} — FAILED ({type(exc).__name__})")
            continue
        total_pages += n
        billed_pages += n if billed else 0
        layer = "textless" if 0 <= chars < 100 else "has-text"
        source = "OCR'd (billed)" if billed else "cached"
        print(f"[{i}/{len(pdfs)}] {pdf.name} — {n} pp, {layer}, {source}")

    print(
        f"\nDone: {len(pdfs) - len(failures)}/{len(pdfs)} files, "
        f"{total_pages} pages ({billed_pages} billed, "
        f"~${billed_pages * COST_PER_PAGE_USD:.2f})"
    )
    for name, err in failures:
        print(f"  FAILED: {name} — {err}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

# Interlude — the literature OCR tool (Mistral) for the literature-researcher

**Status: ⏭ NEXT after compaction (user-ruled 2026-07-22).** Build → use on the
library → compact again → resume the #231 queue (#304 · ch15 · #231 Phase 2 ·
task #10). This file is the whole re-anchor; a fresh session needs only it plus
the referenced RX7 files.

## Why (the incident)

Two literature-researcher dispatches over
`scratch/literature/Sanchez(1982) A review of neutron transport approximations.pdf`
(a textless 56-page scan) died on `API Error: 400 Output blocked by content
filtering policy` — the API's output classifier misfiring on sustained verbatim
transcription of neutron-multiplication literature. The structural cause: with
no text layer, an agent must *generate* the paper's text as model output, and
the output stream is what the filter scores. Text that exists on disk enters as
tool-result *input* and leaves as short selected quotes — a different filter
profile entirely (the whole session carries this content in context without
issue). OCR therefore (1) removes the filter surface, (2) converts per-session
visual-read token cost into a one-time preprocessing step, (3) makes the
library greppable, and (4) emits equations as LaTeX (Mistral document OCR).

## User rulings (2026-07-22, verbatim intent)

- **Mistral for everything** — the whole `scratch/literature/` folder (57
  top-level PDFs), not just the textless subset. "I have no problems with
  Mistral having access to my literature folder. I mean, it IS published work."
- **The scan stays the SSOT.** OCR sidecars are for search, navigation, and
  prose quotation; any **load-bearing equation** is spot-verified against the
  rendered page before it enters a theory page or a solver (the ERR-032-class
  hazard: a plausible-but-wrong reference equation).
- Reuse the existing Mistral work in `~/git/RX7_workshop_manual/`.
- After building AND using the tool, compact again and return to the plan queue.

## The RX7 reuse map (surveyed 2026-07-22)

Canonical tool: `~/git/RX7_workshop_manual/manual_figures/extract.py` (695 ln).
The pieces to lift (line anchors verified):

| Anchor | What | Reuse |
|---|---|---|
| `extract.py:62` | `Mistral(api_key=key)` client init; key auto-loaded from a repo-root `.env` (`MISTRAL_API_KEY=...`) | same pattern; ORPHEUS-root `.env` |
| `:76` | model `mistral-ocr-latest` | same |
| `:85` | `client.ocr.process(**kwargs)` (Document OCR; RX7 passes `include_blocks=True` for figure blocks) | same call; blocks optional for us (no figure-cropping needed) |
| `:91` | `ocr_cached(client, pdf_path, cache_path, pages=None)` — caches the raw response as `<stem>.mocr.json`; **re-runs never re-bill** | lift verbatim-adapted; the idempotency backbone |
| `:454` | `page_label(pj, pno)` — per-page labels for headers | adapt to `## p. <n>` sidecar headers |
| RX7 `CLAUDE.md` | run instructions, `--validate` cheap-first pattern (page 0 only), cost note ~$1/1000 pages | mirror the validate-first CLI flag |

RX7's OpenCV frame-detection / figure-cropping stages are NOT needed — the
literature tool is the slim path: render→OCR→cache→markdown sidecar.

## The tool spec

`tools/ocr_literature.py` (the `tools/` dir exists), run with
`.venv/bin/python` (Host rule). Behavior:

1. **Inputs**: `--input-dir` (default `scratch/literature`), `--glob`,
   `--validate` (OCR page 0 only, print the result, no sidecar), `--force`
   (ignore cache).
2. Per PDF: detect text layer (pypdf extracted-text length over first 3 pages;
   report has-text vs textless — OCR runs on **everything** per the ruling, the
   detection is telemetry only).
3. Call Mistral Document OCR (`mistral-ocr-latest`) via the `ocr_cached`
   pattern → cache `scratch/literature_ocr/<stem>.mocr.json`.
4. Emit `scratch/literature_ocr/<stem>.md`: a small provenance header (source
   PDF path + "OCR sidecar — the SCAN is the SSOT; spot-verify load-bearing
   equations against the rendered page" + journal-page mapping if detectable),
   then per-page sections `## p. <N>` with Mistral's per-page markdown
   (equations arrive as LaTeX). **Never flatten page boundaries** — the corpus
   citation discipline is page-anchored.
5. Print a cost/pages summary at the end.

Housekeeping in the same commit: add `.env` to `.gitignore` (**it is NOT
currently ignored** — verified 2026-07-22); add `mistralai` + `pypdf` (already
present) + `python-dotenv` to the dev environment (document: installed into
`.venv`, not added to `pyproject.toml` production deps — it is tooling).
`scratch/literature_ocr/` is untracked by design (derived artifacts of the
user's untracked library). The KEY: the user copies `MISTRAL_API_KEY` into an
ORPHEUS-root `.env` (a working key exists in the RX7 repo's `.env` — never
print either file).

## Validation anchor (built-in ground truth)

OCR `Sanchez(1982)...pdf` FIRST and diff its pp. 481–482 output against the
visually-verified quotes in
`.claude/plans/phase_i_survey_sanchez_mccormick_1982.md` (the selection-criteria
page: "which computer programs are readily available!", the
optically-large/thin rule, the balance-locality sentences). Agreement ⟹ the
pipeline is trustworthy for prose; then sweep the folder.

## The rules edit (part of the interlude)

`.claude/rules/delegation.md`, the "Briefing a literature pull" section gains:

- Check `scratch/literature_ocr/<stem>.md` FIRST (before any visual page read);
  the sidecar is the search/quote surface, the PDF the verification surface.
- The scan is the SSOT: load-bearing equations are spot-verified against the
  rendered page (Read with `pages=`) before transcription into docs/code.
- Output discipline (the filter lesson): paraphrase + page-cite; keep verbatim
  quotes short; write findings to the deliverable file INCREMENTALLY so a
  mid-run kill preserves progress.

## Sequence for the executing session

1. Build `tools/ocr_literature.py` (lift from the RX7 anchors) + `.gitignore`
   `.env` line + deps into `.venv`.
2. Ask the user to drop `MISTRAL_API_KEY` into the ORPHEUS-root `.env` (blocked
   on the user only here).
3. `--validate` on Sanchez-1982 → full Sanchez run → diff vs the memo.
4. Sweep `scratch/literature/` (57 top-level PDFs; report page count + cost as
   it goes; the cache makes re-runs free).
5. The delegation.md rules edit.
6. Commit (tool + gitignore + rules; sidecars stay untracked). Update this
   plan's status to DONE, tracker task #18 → completed.
7. **Compact; resume the queue**: #304 · ch15 (#276-blocked) · #231 Phase 2 ·
   task #10.

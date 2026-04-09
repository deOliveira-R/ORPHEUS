# QA Test Infrastructure

Last updated: 2026-04-06.

## Running Tests

```bash
.venv/bin/python -m pytest tests/ -v -k "not slow"   # ~500 tests, ~2 min
.venv/bin/python -m pytest tests/ -v                  # +slow (~10 min)
```

## Known Gaps

- **MMS tests**: No manufactured-solution tests for SN, MOC, or diffusion. Biggest L1 gap.
- **L3 validation**: No ICSBEP/IRPhE comparison. Aspirational for educational code.
- **Richardson caching**: Heterogeneous references recomputed every run (~15 min).
- **L0 markers**: No `@pytest.mark.l0` tagging yet.

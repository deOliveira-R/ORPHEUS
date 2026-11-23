"""The tracked files of the concept table, spelled once: the table, the pinned debt and the fingerprint seed."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / "docs" / "architecture"
TABLE_FILE = "docs/architecture/concept_table.toml"
DEBT_FILE = "docs/architecture/concept_debt.toml"
SEED_FILE = "docs/architecture/concept_seed.toml"

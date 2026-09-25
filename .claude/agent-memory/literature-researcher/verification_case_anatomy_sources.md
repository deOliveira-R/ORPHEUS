---
name: verification-case-anatomy-sources
description: V&V benchmark-structure sources extracted for #405 reference redesign (Oberkampf-Trucano 2007, ICSBEP Briggs 2003, Ganapol 2008, ANL-7416, Build systems a la carte, Dolstra) with where they live and which are unobtainable
metadata:
  type: reference
---

Extracted 2026-09-24 for #405 step 2 (reference generation/cache redesign). Memo was in the session scratchpad `refgen_literature.md` (ephemeral); key page cites:

- Oberkampf-Trucano 2007 NEA-6298 invited paper, LOCAL (free oecd-nea.org). Printed p = PDF p + 60. §3.1 four elements (purpose / method-free math description / per-SRQ accuracy+pedigree / documentation incl. arithmetic precision); types 1-4 (analytical, MMS, ODE, PDE); §3.2 comparisons NOT stored in the benchmark; pass = observed order for types 1-2.
- Briggs-Scott-Nouri 2003 NSE 145:1 (from ~/Downloads/NSE zip, now LOCAL): ICSBEP sections 1.0 description / 2.0 evaluation (acceptability, >1% keff unc. unacceptable) / 3.0 spec + benchmark-model keff / 4.0 sample calcs (not peer-reviewed) / App A inputs; revision history, cite edition. Dean 2003 NSE 145:20 also copied, unread.
- Ganapol 2008 NEA 6292, LOCAL (free db-doc2008-1.pdf): "benchmark quality" = 4-5 correct digits, per-benchmark "qualification" table (requested-error ladder), hierarchy exact/near/semi/purely-numerical (PDF p.21), Table A 8-facet classification C1-C6 problem, C7-C8 method (PDF p.23).
- ANL-7416 Supp.2 excerpt BSS-15 LOCAL (scanned, corephysics mirror): 3-level Source Situation -> Problem (ID 15-A1) -> Solution (ID 15-A1-1), each with its own submit/accept dates. §III Guidelines & Format UNREAD.
- Mokhov-Mitchell-Peyton Jones 2020 JFP (LOCAL): verifying / constructive / deep constructive traces; deep = determinism required, no early cutoff; §8.8 self-tracking pessimism for full-language tasks.
- Dolstra 2006 thesis (LOCAL): input-addressed (extensional) vs content-addressed (intensional) store.
- NOT obtainable: Oberkampf-Roy 2010, Knupp-Salari 2003, Roache books (paywalled); SAND2000-1444 free but OSTI resets connections from this host (curl AND WebFetch) — ask user to fetch.

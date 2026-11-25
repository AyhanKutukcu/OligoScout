# Changelog

## 0.5.0 — 2026-11-25

Virtual PCR (scope of Kalendar et al. 2024, virtualPCR; non-human genomes excluded)
- **IUPAC-degenerate primers** in the genome search: `search_terminal_rule` matches a
  template base against the code's set; `oligoscout_panel` routes IUPAC primers in
  anchored mode to the equivalent rule (window = anchor, 0 terminal mismatches).
  Brute-force equal on synthetic data and genome-wide (degenerate RT primers, 646,632
  sites, FM-index = exhaustive scan).
- **`oligoscout_vpcr_scan`** (C++, bit-parallel, exhaustive): any FASTA, many records,
  linear or circular, oligos 8-512 nt, IUPAC in oligo and template, template-N policy,
  3'-terminal rule, optional tolerance of G.T/G.G/G.A mismatches. Identical to the Python
  oracle.
- **`cli/virtual_pcr.py`** and `oligoscout_thermo/vpcr.py`: products in pairs mode or
  all-vs-all fingerprinting (single-primer products included), products across the
  origin and inverse PCR on circular templates, non-binding 5'/3' tails (lower case),
  probes (sites + "inside product"), linked-sequence search, amplicon and template
  sequences, site Tm (multi-mismatch Allawi-Peyret NN) and Ta = min Tm + ln(length).
  Exposed as MCP `oligoscout_virtual_pcr`, the local web API (`vpcr`) and a web form.
- Validation against the paper's Zenodo outputs (`cli/vpcr_validate.py`): six human
  genes 6/6 identical coordinates; degenerate Copia RT primers 11/11 virtualPCR
  amplicons found.

## 0.4.1 — 2026-11-25

LNA
- **LNA duplex thermodynamics** (McTigue, Peterson & Kahn 2004, Table 4, 32NN set):
  every `+N` base adds ΔΔH/ΔΔS for its 5′ and 3′ nearest-neighbour pairs to the
  DNA duplex; salt and concentration are applied as for DNA. New module
  `oligoscout_thermo/lna.py`, table `LNA_MCTIGUE_2004` in `params.py`.
- Validation against the paper's Table 1 (`cli/lna_validate.py`, data in
  `reference_py/validation/`): 98 oligos, RMSD ΔTm 0.886 °C with the measured DNA
  references (paper: 0.88 °C); with OligoScout's DNA model RMSD Tm 1.61 °C, the same
  as the DNA references alone (1.54 °C). LT6a/LT6b excluded: the printed sequence
  and trinucleotide disagree.
- Measured domain is one internal LNA (≥ 2 nt from each end) on a perfect-match
  DNA target. Adjacent, terminal or multiple LNAs and LNA next to a mismatch are
  computed additively and flagged `in_domain: false` with reasons.
- Used by TH-01/TH-03 (`engine.evaluate` accepts `+N` and `/end/` notation),
  PAIR-01, CT-05 (blocker/probe discrimination), MCP `oligoscout_thermo` and
  `oligoscout_check_blocker`, the local web API and the browser probe checker
  (JS port identical to Python, 2000/2000). LNA effects on hairpins, dimers and
  3′-end ΔG are not modelled.

## 0.4.0 — 2026-11-22

Specificity
- **Single-primer products.** `oligoscout_panel` now assembles F→ ←F and R→ ←R
  products in addition to F→ ←R / R→ ←F (`--single-primer on|off`, default on).
  Product rows carry `left_role`, `right_role`, `orientation`. Fixed: R–F mismatch
  masks were rendered with the wrong primer length.
- **Sensitive search mode** (`--mode sensitive`): Primer-BLAST's default terminal
  rule without an exact anchor — every site with ≤ 3 mismatches and ≤ 1 in the
  3′-terminal 5 nt. Equal to a brute-force genome scan (32,546/32,546 ctDNA sites,
  498,315/498,315 Test #81 sites) at 1.1–1.6× the anchored search time.
  Exposed in the web form, the local API (`search_mode`) and MCP.
- **SPC-08 complete:** every GRCh38.p14 sequence outside the 24 primary
  chromosomes (chrM, 260 ALT loci, 254 patches, 166 scaffolds incl. rDNA copies)
  is scanned by the new `oligoscout_contig_scan`; it follows the search mode
  (sensitive = Primer-BLAST terminal rule). QC catalog: 91/91 implemented.

Thermodynamics
- **Dimer end terms:** dangling ends, double dangles and terminal mismatches at
  each helix end, chosen like primer3 ntthal. On 1,024 Test #81 pairs the mean
  hetero-dimer ΔG difference to ntthal fell from +1.80 to −0.01 kcal/mol and the
  −6 kcal/mol flag agreement rose from κ 0.41 to 0.85. Ported identically to
  C++ (74/74 Python↔C++) and the browser engine (2000/2000).
- `salt_crosscheck.py` aligns primer3's C/4 concentration convention: 60/60 Tm
  values identical to oligotm.
- RT-05: fixed the strand mapping of the Sugimoto 1995 DNA/RNA hybrid table
  (now equal to Biopython `R_DNA_NN1`).

Integration and packaging
- MCP server (15 tools, stdio + Streamable HTTP + REST/OpenAPI), `search_mode`
  on design and specificity tools, OAuth 2.1 for public clients.
- Dockerfile, root LICENSE, VERSION, YAML QC catalog, ready-made primer catalog.

Benchmarks
- `OligoScout/benchmarks/scorecard/`: thermodynamics vs primer3, brute-force
  blind-spot analysis, sensitive-mode verification, repeated timings on native
  Linux, consolidated `SCORECARD.md`.
- NCBI Primer-BLAST web, 31 pairs: all 654 products that satisfy Primer-BLAST's
  own rule are found by the sensitive mode (653 on the primary chromosomes, 1 on
  a chr22 unlocalized scaffold via SPC-08); 87 more that Primer-BLAST missed.
- PrimerBank, 121 experimentally tested mouse pairs: TH-03 at 60 °C annealing
  fails 43/101 sequence-confirmed pairs (0/101 at 58 °C) — uncalibrated
  threshold, to be revisited with wet-lab data (REP-03).

Known gaps
- Non-human genomes and wet-lab calibration are out of scope.

## 0.3.0
- Thermodynamic reference engine, 90/91 QC rules, Test #81 benchmark, web UI.

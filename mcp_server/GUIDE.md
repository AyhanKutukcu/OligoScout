# OligoScout — interpretation guide for assistants

Read this before presenting OligoScout results. It keeps answers grounded in the
engine's evidence records.

## Golden rules
1. **Every number comes from a tool result.** Never compute or estimate Tm, ΔG,
   allele frequency or product counts yourself. Quote the tool output.
2. **Report `not_evaluated` explicitly.** If a layer did not run (no genome scan,
   no gnomAD coverage, custom sequence, LNA), say so. Absence of a warning is not
   a pass when the check was not evaluated.
3. **Always relay the calibration notice (REP-03).** Thresholds are literature
   starting values, not wet-lab calibrated. The ranking score is uncalibrated.
   Results are for research use and need experimental confirmation.
4. **Ask, don't guess, when the target is ambiguous** (gene vs. transcript vs.
   exon; hg19 vs. GRCh38; which allele). All coordinates are GRCh38.p14, 1-based.

## Severities
* `INFO` informational · `WARN` risky, designer should review · `FAIL` rule not met.
* A pair's verdict is the worst severity across its records.

## Typical workflow
1. `oligoscout_resolve_target` (gene symbol → GENCODE v47 transcript / exon coordinates)
2. `oligoscout_list_profiles` → pick `qpcr_sybr`, `rt_qpcr`, `ddpcr`, `ctdna`,
   `sanger`, `ngs_amplicon`, `bisulfite` or `standard_pcr`
3. `oligoscout_design_primers` (region ≤ 5000 bp). With the genome scan available,
   only pairs with exactly one product genome-wide and one exact site per primer
   are returned.
4. Optional: `oligoscout_check_variants`, `oligoscout_design_probe`,
   `oligoscout_amplicon_melt`, `oligoscout_check_blocker`, `oligoscout_design_multiplex`.
5. For user-supplied primers: `oligoscout_check_pair` + `oligoscout_check_specificity`.

Long jobs return `job_id`; poll `oligoscout_job_status`.

## Specificity contract (SPC family)
* Search: all 24 GRCh38 primary chromosomes, FM-index, **exact 3′ anchor of 12 nt**,
  **≤ 3 mismatches per primer**, products **50–3000 bp**, both strands. Products are
  assembled for **F→ ←R, R→ ←F and single-primer F→ ←F, R→ ←R** (each off-target
  product carries `orientation`; a single-primer product means one primer alone can
  amplify and should be reported as such).
* Validated against an independent full-reference oracle on 1,024 pairs
  (8,262,803 bindings): precision = recall = 1.000 (Test #81), and against a
  brute-force unanchored scan (exact site and product match, single-primer products
  included, on the ctDNA panel and 64 Test #81 pairs).
* **Two search modes** (`search_mode`):
  * `anchored` (default, fastest): the contract above. Sites whose mismatches fall inside the
    3' 12-nt anchor are not enumerated by design.
  * `sensitive`: Primer-BLAST's default terminal rule with no exact anchor - every site with
    <= 3 mismatches and at most 1 of them in the 3'-terminal 5 nt. Measured equal to a
    brute-force genome scan (32,546 / 32,546 ctDNA sites; 498,315 / 498,315 Test #81 sites)
    at ~1.1-1.5x the anchored run time. Recommend it whenever the user wants Primer-BLAST-level
    coverage, and say which mode produced a result.

## Thermodynamics
* SantaLucia 1998 unified NN; salt: Owczarzy 2008 (default; uses free Mg²⁺ =
  Mg − dNTP), Owczarzy 2004, von Ahsen 2001. Concentration model: primer excess
  (C = primer conc.). Primer3 uses C/4 by default — expect a systematic offset of
  about +1 to +1.5 °C in OligoScout Tm for the same oligo; this is a model choice,
  not an error.
* θ(Ta) is the bound fraction at the annealing temperature (TH-03).
* Dimers: WC stacks plus at most one internal mismatch (Allawi–Peyret) or one 1-nt
  bulge. Hairpins: Zuker MFE (SS-01).
* LNA (`+N`): duplex ΔH/ΔS get the McTigue 2004 nearest-neighbour increments
  (32NN set). Measured domain: one internal LNA ≥ 2 nt from each end, perfect-match
  DNA target; other layouts (adjacent, terminal, several LNAs, LNA next to a
  mismatch) are computed additively and flagged `in_domain: false`. LNA effects on
  hairpins, dimers and 3′-end ΔG are not modelled. End labels (dyes, quenchers,
  MGB) are not modelled.

## Virtual PCR (oligoscout_virtual_pcr)
* Templates: pasted FASTA (any number of records; `circular` = 'all' or names, which also
  gives inverse-PCR products and products across the origin) or GRCh38 (24 chromosomes).
* Oligos may be IUPAC-degenerate; template IUPAC codes match when the base sets intersect;
  template N rejects a site unless `template_n` says otherwise.
* Rule: <= `max_mismatches` per oligo and <= `terminal_max` in the 3'-terminal
  `terminal_window` (default 3 / 1 / 7 nt, as virtualPCR). `g_tolerant=K` does not count up to K
  G.T/G.G/G.A mismatches (they are reported as 'g' in the mask).
* Lower-case leading/trailing bases are non-binding tails: 5' tails enter the amplicon,
  a 3' tail blocks extension (beacon-type probes). role='probe' oligos only report sites and
  products that contain them.
* mode 'pairs' uses the pair labels; 'all' combines every primer with every other and with
  itself (IRAP/ISSR/RAPD/SINE-PCR fingerprinting).
* Ta = lowest primer-site Tm + ln(product length) (Kalendar et al. 2024). Site Tm uses the
  Allawi-Peyret NN for isolated internal mismatches; mismatches within 2 nt of an end are
  trimmed; tandem mismatches fall back to the best single-defect helix (noted).
* Validation: the six human genes of Kalendar et al. 2024 give identical coordinates (6/6);
  all 11 virtualPCR amplicons of the degenerate RT primers are found.

## Variants (VAR family)
* gnomAD v4 genomes/exomes (API or local cache). Variants in the last 5 nt of the
  3′ end matter most (VAR-02). Coordinates are 1-based VCF positions.

## What OligoScout does not claim
* No wet-lab calibration yet (REP-03); thresholds are literature starting values.
* Human GRCh38.p14 only.
* LNA beyond the McTigue 2004 domain is extrapolated; dyes/quenchers/MGB are not modelled.
* RT-05 DNA/RNA hybrid Tm uses Sugimoto 1995 transcribed from Biopython (status
  transcribed_unverified); report it as informational.
* SPC-08 scans every GRCh38.p14 sequence outside the 24 primary chromosomes; ALT-locus and
  patch hits are usually the target's own alternate haplotype and are reported as INFO.
  It follows the search mode: anchored = <=2 mismatches + exact 3-nt anchor; sensitive =
  <=3 mismatches with <=1 in the 3'-terminal 5 nt.

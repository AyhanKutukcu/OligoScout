# OligoScout validation & benchmark scorecard

_generated 2026-11-24T12:55:28Z_

## 1. Specificity search vs. independent oracle (Test #81)

| Method | Binding P | Binding R | Product R | Search s | Max RSS GiB |
|---|---:|---:|---:|---:|---:|
| Independent full-scan oracle | 1.000000 | 1.000000 | 1.000000 | 194.1 | 0.24 |
| OligoScout | 1.000000 | 1.000000 | 1.000000 | 111.8 | 0.65 |
| Bowtie 1 + contract filter | 1.000000 | 1.000000 | 1.000000 | 148.6 | 3.65 |
| BWA-aln + contract filter | 1.000000 | 1.000000 | 1.000000 | 141.1 | 5.03 |
| BLASTN-short e=1100 + contract filter | 1.000000 | 0.999175 | 0.999995 | 3065.4 | 10.98 |

## 2. Anchor blind spot (exhaustive unanchored <=3-mismatch scan)

**ctdna_panel** (20 pairs, 22 primers)

- sites: 39,049 with <=3 mm; anchored (OligoScout contract) 4,931; Primer-BLAST-rule relevant 32,546; blind spot 27,615
- products: OligoScout 26; Primer-BLAST rule 728; extra with <=2 mm/site 38; with <=1 mm/site 0; single-primer products (F-F/R-R) found by oligoscout_panel 1 of 1 in the contract
- pairs with any extra product: 13/20 (<=2 mm/site: 3/20)

**test81_first64_balanced** (64 pairs, 128 primers)

- sites: 512,082 with <=3 mm; anchored (OligoScout contract) 12,141; Primer-BLAST-rule relevant 498,315; blind spot 486,174
- products: OligoScout 1,683; Primer-BLAST rule 95,153; extra with <=2 mm/site 20,177; with <=1 mm/site 4,776; single-primer products (F-F/R-R) found by oligoscout_panel 1 of 1 in the contract
- pairs with any extra product: 15/64 (<=2 mm/site: 11/64)

## 2b. Sensitive mode (Primer-BLAST terminal rule) vs. brute force

| Panel | sites brute force | sites OligoScout | products brute force | products OligoScout | exact |
|---|---:|---:|---:|---:|---|
| ctdna | 32,546 | 32,546 | 728 | 728 | yes |
| test81_64 | 498,315 | 498,315 | 95,153 | 95,153 | yes |

## 2c. NCBI Primer-BLAST web vs. OligoScout sensitive mode

- 31 pairs (20 ctDNA test designs (DESIGN_000-019) + 11 ready-made catalog pairs); Primer-BLAST listed 2,229 products, 654 satisfy the terminal rule on both sites (88 listed products violate it)
- rule products also found by OligoScout: 653 (primary chromosomes 653/653 (100%); with SPC-08 extra contigs 654/654 (100%))
- OligoScout-only rule products (Primer-BLAST seeding misses): 87

## 2d. Repeated timings (native Linux VM)

_AMD Ryzen 9 9950X (VM, 2 vCPU, 3 GiB RAM), Ubuntu 22.04, native Linux kernel; PPFM index read from the connected Windows folder mount; single search thread_

| Run | n | search s (mean ± SD) | max RSS MiB | products |
|---|---:|---:|---:|---:|
| ctdna_anchored | 5 | 26.16 ± 0.12 | 660.5 | 27 |
| ctdna_sensitive | 5 | 28.98 ± 0.12 | 660.6 | 728 |
| t81_64_anchored | 5 | 26.69 ± 0.15 | 660.6 | 1,684 |
| t81_64_sensitive | 5 | 42.96 ± 0.08 | 662.5 | 95,153 |
| test81_1024_anchored | 4 | 33.73 ± 0.41 | 684.5 | 425,055 |

## 3. Thermodynamics vs. primer3 (oligotm / ntthal), 2,048 primers

| Comparison | mean diff | SD | max abs | Pearson r |
|---|---:|---:|---:|---:|
| A_tm_1M_aligned_conc | +0.000 | 0.000 | 0.000 | 1.0 |
| C_tm_same_salt_model_aligned_conc | -0.000 | 0.000 | 0.000 | 1.0 |
| B_tm_default_each_tool | +0.995 | 0.228 | 1.833 | 0.9052 |
| D_heterodimer_dG37 | -0.007 | 0.462 | 3.595 | 0.9617 |
| E_selfdimer_dG37 | +0.089 | 0.709 | 4.490 | 0.9446 |
| F_hairpin_dG37 | +0.265 | 0.398 | 2.348 | 0.3453 |

Dimer model: nn_single_internal_mismatch_or_bulge+end_terms (dangling ends, double dangle, terminal mismatch). Before end terms (nn_single_internal_mismatch_or_bulge (no end terms)): hetero-dimer mean diff +1.80, kappa 0.4094; now kappa 0.851 (self-dimer 0.8207).

## 3b. PrimerBank experimentally tested mouse pairs

- 121 PrimerBank mouse pairs with validation pages (89 genes): 101 sequence-confirmed, 9 expected transcript not the top BLAST hit, 11 without a BLAST result
- OligoScout worst verdict on sequence-confirmed pairs (rt_qpcr, Ta 60 C): {'WARN': 56, 'FAIL': 45}; FAIL mostly TH-03 (43 pairs)
- TH-03 failures on working pairs by annealing temperature: 60 C: 43/101, 58 C: 0/101, 56 C: 0/101
- Tm OligoScout - PrimerBank: mean +1.59 C (SD 0.33)
- PrimerBank recommends 60 C annealing. At 60 C OligoScout's TH-03 (bound fraction < 0.80) fails pairs that were sequence-confirmed, because OligoScout's Tm under 50 mM Na+/1.5 mM Mg2+ sits near 60-64 C for these primers; TH-02 (salt-model spread > 2 C) warns on every primer. Both thresholds are uncalibrated (REP-03) and this sample argues for relaxing them.

## 4. Software integrity

- unit tests: 279 OK; C++ single-primer + assembler equivalence tests PASS
- end-to-end regression: {'pass': 20, 'skip': 1, 'fail': 0, 'evidence': 'reference_py/out/regresyon/REGRESYON_RAPORU.md'}
- Python vs C++: 74 sequences, 0 differences
- parameter tables vs primer3: 371/373 identical, 0 flagged deviations
- salt layer vs oligotm: 60/60 identical
- QC catalog: 91 rules, status {'etiketli_uygulandi': 91}

## 5. Not yet benchmarked

- Wet-lab validation (qPCR efficiency, specificity by melt/gel/NGS) - required before any calibrated claim
- Head-to-head design quality vs Primer3/Primer-BLAST/PrimerQuest outputs on the same targets with a wet-lab outcome (PrimerBank gives outcomes only for its own mouse designs; QC false-alarm rate on them is in section 3b)
- Multiplex panel performance vs PrimerPlex/Ultiplex/olivar on wet-lab data
- LNA-modified oligos (needs the McTigue 2004 LNA NN table) and non-human genomes (index + annotations)

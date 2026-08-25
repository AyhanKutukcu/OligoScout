# Test #81 publication benchmark — 1024 primer pairs

Contract: 24-nt primers; exact 12-nt 3′ anchor; at most 3 mismatches per primer; 50–3000 bp products; GRCh38.p14 primary assembly.

| Method | Binding P/R/F1 | Product P/R/F1 | Intended recall | Search wall | End-to-end wall | Max stage RSS |
|---|---:|---:|---:|---:|---:|---:|
| Independent full-scan oracle | 1.000000 / 1.000000 / 1.000000 | 1.000000 / 1.000000 / 1.000000 | 1.000000 | 194.11 s | 379.48 s | 2,963,900 KiB |
| OligoScout | 1.000000 / 1.000000 / 1.000000 | 1.000000 / 1.000000 / 1.000000 | 1.000000 | 111.80 s | 255.91 s | 2,963,684 KiB |
| Bowtie 1 + contract filter | 1.000000 / 1.000000 / 1.000000 | 1.000000 / 1.000000 / 1.000000 | 1.000000 | 148.64 s | 479.53 s | 3,832,160 KiB |
| BWA-aln + contract filter | 1.000000 / 1.000000 / 1.000000 | 1.000000 / 1.000000 / 1.000000 | 1.000000 | 141.14 s | 419.17 s | 5,271,468 KiB |
| BLASTN-short e=1100 + contract filter | 1.000000 / 0.999175 / 0.999588 | 1.000000 / 0.999995 / 0.999998 | 1.000000 | 3065.39 s | 3558.70 s | 11,518,096 KiB |

Bowtie 1, BWA-aln, and BLASTN-short are general alignment/search baselines. Their raw hits were independently revalidated against the reference and filtered to the same primer-binding contract. They are not described as primer-native tools.

The independent oracle performs a direct full-reference scan. Product-level uncertainty is reported with 10,000 fixed-seed bootstrap resamples using the primer pair—not individual hits—as the resampling unit.

Timing caveat: large outputs and temporary sort files reside on the Windows E: drive through WSL DrvFS. These wall times describe this machine and storage path; accuracy/set metrics are unaffected. Search wall time is reported separately from an end-to-end sum of all sequential measured stages. BLAST search wall time is the sum of sequential query chunks and its RSS is the maximum observed chunk RSS; adaptive smaller chunks may be used for resource-heavy queries and are retained in the audit trail. OligoScout's native search stage also materializes binding and product tables; external tools first materialize raw alignments and then undergo contract normalization and product assembly. Interpret both the search-phase and end-to-end columns with these scope differences in view.

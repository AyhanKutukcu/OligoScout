# Test #81: publication-scale external benchmark

This directory contains the compact, reviewable artifacts for OligoScout Test
#81. The frozen run identifier is `test81_publication_20260825_031353`; the
benchmark used commit `c89df254df1a9fdda6bb5e16b98166cf1b43425c` as its
source baseline.

## Experimental contract

- Reference: normalized GRCh38.p14 primary assembly, 24 sequences
- Reference SHA-256: `d5fdff32d78baf7eb978b0edd2e9e28e065a691f7797cf53d8278e621ba21f5e`
- Panel: 1,024 Primer3-designed pairs / 2,048 primers
- Primer length: 24 nt
- Exact 3-prime anchor: 12 nt
- Mismatches: at most 3 per primer
- PCR product interval: 50–3000 bp
- Search threads: 1
- Difficulty design: four balanced anchor-frequency strata of 256 pairs
- Statistical unit: primer pair
- Bootstrap: 10,000 fixed-seed resamples

The first 16, 64, 128, and 256 pairs are balanced nested subsets because the
four strata were interleaved deterministically. Intended loci in the full panel
span all 24 primary chromosomes.

## Independent truth set

The oracle does not consume OligoScout, Bowtie, BWA, or BLAST hits. It scans the
reference directly, uses the exact 12-nt anchor only to enumerate candidates,
and independently verifies the complete primer sequence and mismatch count.
Products are then assembled from same-chromosome, inward-facing binding sites.

The oracle produced:

- 8,262,803 binding sites;
- 425,055 products;
- 1,024/1,024 intended products.

## Results

| Method | Binding precision | Binding recall | Product precision | Product recall | Intended recall | Search (s) | End-to-end (s) | Search max RSS (KiB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OligoScout | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 111.80 | 255.91 | 680,804 |
| Bowtie 1 + filter | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 148.64 | 479.53 | 3,832,160 |
| BWA-aln + filter | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 141.14 | 419.17 | 5,271,468 |
| BLASTN-short + filter | 1.000000 | 0.999175 | 1.000000 | 0.999995 | 1.000000 | 3065.39 | 3558.70 | 11,518,096 |

OligoScout, Bowtie 1, and BWA-aln were exactly equal to the oracle at both
binding and product level. BLASTN-short omitted 6,814 valid bindings and two
non-intended products, with no additional contract-valid predictions. The two
missing products belonged to `Q4_HIGH` pairs 851 and 955. All methods recovered
all 1,024 intended products.

BLAST was executed as 31 original 64-query chunks followed by eight adaptive
8-query chunks for the final resource-heavy group. The original 64-query FASTA
and its lossless eight-way split are retained under `provenance/`; concatenating
the split files is byte-identical to the original. Chunking changed the resource
schedule, not the query set or search contract.

## Timing interpretation

Search and end-to-end times are reported separately. OligoScout's native search
stage materializes binding and product tables. General-purpose tools first
materialize raw alignment/search output and then undergo reference validation,
contract filtering, sorting, and product assembly.

The run used a 12th Gen Intel Core i5-12500H under WSL2. Large outputs and sort
temporary files were written through the Windows E: DrvFS mount. Index
construction was excluded. These single-run numbers characterize this frozen
configuration only. Publication-level general performance claims should add
repeated runs with an explicit cache policy, preferably on a native Linux file
system and an independent machine.

## Directory contents

```text
benchmarks/test81/
├── README.md
├── inputs/
│   ├── test81_panel_1024.tsv
│   └── test81_panel_1024.metadata.json
├── results/
│   ├── FINAL_SUMMARY.tsv
│   ├── test81_results.tsv
│   ├── test81_report.md
│   ├── test81_product_f1.svg
│   ├── test81_search_runtime.svg
│   └── comparisons/<method>/
├── provenance/
│   ├── environment.public.tsv
│   └── adaptive_chunking/
├── src/
├── tests/
├── test81_*.py
├── *.sh
└── SHA256SUMS
```

The published `SHA256SUMS` manifest covers every compact artifact in this
directory except the manifest itself. Multi-gigabyte raw SAM/BLAST files,
complete binding/product tables, reference FASTA files, and generated indexes
are intentionally not committed.

## Reproduction scope

The archived drivers preserve the exact executed workflow. They expect a local
GRCh38 primary FASTA, PPFM shards, and the pinned external-tool indexes described
by `test81_common.sh`; these large assets are excluded from Git. Override
`TEST81_REPO` and `TEST81_OUTPUT_ROOT` when using different local locations.

From a prepared environment, run the 128-pair gate first:

```bash
bash benchmarks/test81/setup_e_drive_and_pilot.sh
```

After `TEST81_PILOT_COMPLETE YES`, launch the full panel with the run tag printed
by the pilot:

```bash
TEST81_RUN_TAG=<run-tag> \
    bash benchmarks/test81/launch_full_1024.sh
```

The tiny fixtures under `tests/` exercise plus/minus coordinate reconstruction,
SAM and BLAST normalization, independent oracle enumeration, and product
assembly without requiring GRCh38.

# OligoScout

OligoScout is research-oriented C++20 software for exact and approximate
matching of short DNA patterns using persistent bidirectional FM-index
structures. Its principal biological workflow is strand-aware PCR
primer-pair analysis, including binding-site discovery, off-target product
assembly, multiplex evaluation, mismatch profiling, and reproducible
reference-scale validation.

> **Project status:** pre-1.0 research software. The core libraries,
> validation tools, and benchmark drivers are implemented and tested. The
> top-level CLI currently exposes the FM-index demonstration command. Its
> `build`, `search`, and `search-genome` commands are reserved interfaces and
> are not implemented yet.

## Scope

OligoScout is designed for:

- exact short-DNA pattern search;
- bounded-mismatch and anchor-constrained search;
- bidirectional FM-index traversal;
- strand-aware genomic coordinate reporting;
- PCR primer-pair and off-target product assembly;
- multiplex primer cross-product analysis;
- biological mismatch and thermodynamic annotation;
- structured TSV, JSON, and HTML reporting;
- persistent chromosome-sharded PPFM indexes.

OligoScout is **not currently** a general-purpose sequencing-read aligner.
It does not claim support for arbitrary gapped alignment, CIGAR/MAPQ
generation, soft clipping, splice-aware RNA alignment, or multiple-sequence
alignment.

## Architecture

```text
Reference FASTA
      |
      v
Packed BWT / FM-index / bidirectional FM-index
      |
      v
Persistent chromosome-sharded PPFM indexes
      |
      v
Exact, sensitive, candidate and adaptive short-pattern search
      |
      v
Strand-aware binding sites
      |
      v
Primer-pair and multiplex product assembly
      |
      v
Biological features, thermodynamics and risk annotation
      |
      v
TSV / JSON / HTML reports
```

## Requirements

- Linux or WSL2;
- a C++20 compiler;
- CMake 3.25 or newer;
- Git, including submodule support;
- GNU Make for the pinned Primer3 dependency.

Optional benchmark and data-preparation workflows require additional tools.
Generated indexes, whole-genome references, benchmark outputs, local tool
installations, and virtual environments are intentionally excluded from this
repository.

## Clone and build

```bash
git clone --recurse-submodules \
    https://github.com/AyhanKutukcu/OligoScout.git

cd OligoScout

make -C third_party/primer3/src -j4

cmake --preset release

cmake \
    --build \
    --preset release \
    --target primerpair-bifm \
    -j4
```

The internal CMake target remains `primerpair-bifm` for source compatibility.
The generated user-facing executable is:

```text
build/release/oligoscout
```

## Current CLI

```bash
./build/release/oligoscout --help
./build/release/oligoscout demo
```

The following names are visible in help output but are placeholders at this
stage and return a clear “not implemented” error without reading a FASTA file
or creating an index:

- `build`
- `search`
- `search-genome`

The implemented research workflows are presently exercised through the C++
libraries, test executables, dedicated tools, and benchmark drivers.

## Tests

```bash
ctest \
    --test-dir build/release \
    --output-on-failure
```

The publication baseline contains 78 Release regression tests covering the
FM-index layers, persistent indexes, search modes, chromosome sharding,
primer-pair assembly, multiplex behavior, biological annotation, reports,
and supporting data structures.

## External benchmark baseline

The frozen Test #79 comparison used 32 primer pairs (64 primers), GRCh38.p14
primary chromosomes, a 50–3000 bp product interval, at most three mismatches
per primer, and an exact 12-nt 3-prime anchor for the matched contract.

| Result | OligoScout / historical PrimerPair-Search | BLAST+ matched contract |
|---|---:|---:|
| Products | 65,467 | 65,454 |
| Known-target recall | 32/32 | 32/32 |
| Shared products | 65,454 | 65,454 |
| Method-only products | 13 | 0 |
| Matched-contract Jaccard | 0.999801427 | 0.999801427 |

The 13 OligoScout-only products were independently revalidated. They
represented 10 unique BLAST binding sites omitted by the original BLAST save
threshold: no sites were recovered with `evalue 1000`, all 10 were recovered
with `evalue 1100`, and the observed HSP E-value was approximately 1060.

In that single-thread run, OligoScout completed in 17.02 seconds with a
maximum RSS of 680,172 KB, while BLAST+ completed in 179.37 seconds with a
maximum RSS of 11,412,580 KB. These values describe that specific benchmark
configuration only and are not a universal speedup or memory claim.

MFEprimer was recorded as environment-limited: its index build did not
complete and its search benchmark was not executed. No MFEprimer performance
or accuracy claim is made. Later multi-aligner experiments include
resource-limited profiles and are not presented here as completed comparative
claims.

## Reproducibility and repository contents

The repository tracks source code, headers, tests, benchmark drivers, small
test inputs, configuration, and compact trained model artifacts. It does not
track:

- GRCh38 FASTA files or source archives;
- generated PPFM, Bowtie, BWA, HISAT2, minimap2, or Jellyfish indexes;
- benchmark result trees, logs, attestations, or backups;
- build trees and Python virtual environments;
- locally installed third-party binaries.

Primer3 is recorded as a Git submodule at commit
`345cba9ca22fdb8f9f90da9b57948b81b4866330`.

## Naming note

The current user-facing project and executable name is **OligoScout**.
Historical benchmark artifacts and some internal identifiers retain the
earlier `PrimerPair-Search`, `PrimerPair-BiFM`, or `primerpair` names to
preserve reproducibility and avoid an algorithm-changing namespace migration.

## License

This project is distributed under the MIT License. See [LICENSE](LICENSE).

## Maintainer

Ayhan Kutukcu

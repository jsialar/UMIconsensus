# minihap_workflow

A [Nextflow](https://www.nextflow.io/) pipeline for generating **UMI-based consensus sequences (minihaplotypes)** from Oxford Nanopore (ONT) amplicon sequencing data.

The pipeline takes raw ONT reads, corrects errors by collapsing reads that share the same Unique Molecular Identifier (UMI) into a single high-accuracy consensus sequence, and calls the dominant haplotype(s) per amplicon target.

> Manifest: `jsialar/UMIconsensus` · version `v0.9.0` · author Jaren Sia

## Overview

```
fastq (per barcode)
      │
      ▼
  MERGE_FASTQ ──► TRIM_FASTQ ──► MAP_READS ──► GENERATE_CONSENSUS
   (cat)          (cutadapt)     (minimap2)          │
                                                     ├─ PARSE_BED       split bed into per-target regions
                                                     ├─ GROUP_UMI       group reads by UMI (umi_tools)
                                                     ├─ SPLIT_GROUP     split BAM per UMI family (samtools)
                                                     ├─ CALL_CONSENSUS  build consensus per family
                                                     └─ PROCESS_CONSENSUS  filter noise → haplotype + quality tables
```

### Pipeline steps

| Step | Module | Tool | Description |
|------|--------|------|-------------|
| Merge | `MERGE_FASTQ` | `cat` | Concatenate all `*.fastq.gz` files within each `barcode*` folder into one file per sample. |
| Trim | `TRIM_FASTQ` | `cutadapt` | Trim sequencing adapters, filter by read length, and move the UMI into the read name. |
| Map | `MAP_READS` | `minimap2` + `samtools` | Align trimmed reads to the reference genome and produce a sorted, indexed BAM. |
| Group UMI | `GROUP_UMI` | `umi_tools` (modified) | Group reads into UMI families per amplicon target, filtered by family size. |
| Split | `SPLIT_GROUP` | `samtools split` | Split the grouped BAM into one BAM per UMI family. |
| Call consensus | `CALL_CONSENSUS` | Python + `samtools consensus` | Call a consensus base at each SNP position for every UMI family, giving one minihaplotype per family. |
| Process consensus | `PROCESS_CONSENSUS` | Python | Drop low-frequency (noise) consensus calls and emit the dominant haplotypes with per-base quality. |

## Requirements

- [Nextflow](https://www.nextflow.io/) `>=23.04.2`
- [Docker](https://www.docker.com/) (recommended), or the tools below installed locally:
  - `minimap2` 2.30, `samtools` 1.21
  - Python with `pysam` 0.22.1, `pandas` 2.2.2, `cutadapt` 5.1, `umi_tools` 1.1.6, `dnaio`

A prebuilt container image (`m113track/minihap`) is used by the `standard` profile. To build it yourself:

```bash
docker build -t m113track/minihap Docker/
```

## Input

Point `--input` at the **parent folder** containing per-barcode subfolders of ONT reads:

```
<input>/
├── barcode01/
│   ├── *.fastq.gz
│   └── ...
├── barcode03/
│   └── *.fastq.gz
└── ...
```

You also need:
- `--bed` — BED file of amplicon target sites (tab-separated: `chrom  start  end  name`).
- `--snplist` — CSV listing the SNPs within the amplicon sites (`Chromosome,Start,End,Name`).
- `--reference` — genome reference, indexed (`*.mmi`) or FASTA.

Example files are provided under [data/](data/) and [tests/data/](tests/data/).

## Usage

```bash
nextflow run main.nf \
  --input   /path/to/fastq_pass \
  --bed     data/0.01Mslack_refined_cleaned_filtered.bed \
  --snplist data/nonpolysnp_nonindel.csv \
  --reference /path/to/reference.mmi \
  --output  results/myrun \
  -profile  standard
```

### Profiles

- `standard` — runs each process inside the `m113track/minihap` Docker container (default containerised setup).
- `local` — Docker disabled; expects the tools to be available on `PATH`.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *(required)* | Parent folder containing `barcode*/` subfolders of `*.fastq.gz`. |
| `--output` | `null` | Output folder. |
| `--reference` | `.../GRCh38_no_alt_analysis_set.mmi` | Reference genome (`*.mmi` or `*.fasta`). |
| `--bed` | *(required)* | BED file of amplicon target sites. |
| `--snplist` | *(required)* | CSV of SNPs within the amplicon sites. |
| `--minlength` | `800` | Minimum read length (excluding adapter). |
| `--maxlength` | `1200` | Maximum read length (excluding adapter). |
| `--save_untrimmed` | `false` | Keep reads without adapter sequences in a separate file. |
| `--method` | `adjacency` | UMI correction method: `adjacency`, `directional`, or `cluster`. |
| `--familysizethreshold` | `5` | Minimum UMI family size to keep. |
| `--keep_nonconsensus` | `false` | Keep consensus sequences containing non-substitution calls (`N`/`D`/`I`/`U`) — see [Consensus base codes](#consensus-base-codes). |
| `--noise_threshold` | `0.05` | Drop haplotypes below this fraction of the most abundant one. |
| `--minimap2_param` | `-ax map-ont -k 13 --MD` | minimap2 alignment parameters. |
| `--additional_outputs` | `false` | Also publish per-family consensus FASTQs and per-family BAMs. Only honoured by the [generate_consensus_onwards.nf](generate_consensus_onwards.nf) entry point; `main.nf` ignores it. |
| `--debug` | `false` | Keep the Nextflow work directory (skip cleanup). |

See [nextflow_schema.json](nextflow_schema.json) for the full schema.

## Output

Results are published under `<output>/<barcode>/<module>/`:

```
<output>/<barcode>/
├── trim/
│   ├── <barcode>.json                  cutadapt report (adapters found, length filtering)
│   └── <barcode>.untrimmed.fastq.gz    reads with no adapter — only with --save_untrimmed
├── aligned/
│   ├── <barcode>.bam                   minimap2 alignments, coordinate sorted
│   └── <barcode>.bam.bai
├── group_umi/
│   └── <barcode>_umisizes.tsv          UMI family size distribution, all targets concatenated
└── call_consensus/
    ├── <barcode>_consensustable.csv    one row per distinct minihaplotype per target
    ├── <barcode>_haps.csv              dominant haplotypes surviving the noise filter
    └── <barcode>_quals.csv             per-base quality for those dominant haplotypes
```

File formats (all CSVs are written without a header row):

| File | Columns | Notes |
|------|---------|-------|
| `<barcode>_umisizes.tsv` | UMI family size distribution per target | Concatenated across all amplicon targets. |
| `<barcode>_consensustable.csv` | `minihaplotype, umifamilysize, count, target` | `count` is how many UMI families gave that exact `(minihaplotype, umifamilysize)` pair, so one haplotype spans several rows when its families differ in size. Concatenated across targets. |
| `<barcode>_haps.csv` | `<target>_<idx>, minihaplotype` | Haplotypes at or above `--noise_threshold` of the most abundant one for that target, numbered from 0. |
| `<barcode>_quals.csv` | `<target>_<idx>, q1, q2, …` | One Phred score per SNP position, taking the **highest** score seen at that position across all UMI families supporting the haplotype. Row *idx* matches the same *idx* in `_haps.csv`. |

`CALL_CONSENSUS` also writes a per-target `<target>.csv` and `<barcode>.<target>.fastq.gz` inside the work directory. The CSVs are concatenated into `<barcode>_consensustable.csv` above; the per-family FASTQs are **not** published by `main.nf` (see `--additional_outputs` below).

### Running from an existing BAM

[generate_consensus_onwards.nf](generate_consensus_onwards.nf) is a second entry point that skips merge/trim/map and starts from an already-aligned BAM:

```bash
nextflow run generate_consensus_onwards.nf \
  --bam '/path/to/<barcode>.bam' \
  --bed data/0.01Mslack_refined_cleaned_filtered.bed \
  --snplist data/nonpolysnp_nonindel.csv \
  --output results/myrun
```

It publishes the same `group_umi/` and `call_consensus/` outputs, and it is the **only** entry point that honours `--additional_outputs`, which adds:

```
<output>/<barcode>/
├── fastq_consensus/   <barcode>.<target>.fastq.gz   one record per UMI family, with per-SNP qualities
└── splitted_bam/      P*.bam + P*.bam.bai           one BAM per UMI family
```

### Consensus base codes

A consensus "sequence" is a **minihaplotype**, not a contiguous stretch of the genome: it has exactly one character per SNP listed in `--snplist` for that amplicon, in the order those rows appear in the CSV. Character *i* is the call for SNP *i*, which is what makes the strings comparable across UMI families.

Most characters are ordinary bases, but four codes mark positions that are not a clean substitution:

| Code | Meaning | How it arises |
|------|---------|---------------|
| `A` `C` `G` `T` | Substitution call | `samtools consensus` returned a single base for the position. |
| `N` | Ambiguous | The position is covered, but the reads disagree too much for a confident call. Emitted by `samtools consensus` itself. |
| `D` | Deletion | The position is covered but the reference base is deleted in this UMI family (`samtools consensus --show-del yes` returns `*`). |
| `I` | Insertion | `samtools consensus` returned more than one base for the single-base window, meaning an insertion immediately **after** the position. Insertions are anchored to the preceding reference base, so an insertion sitting just *before* the SNP is attributed to the previous position and does not raise `I` here. |
| `U` | Unknown / no coverage | This UMI family has no reads spanning the position at all — truncated or soft-clipped reads, or a family with no alignments. `samtools consensus` emits nothing, as distinct from the covered-but-deleted `D` case. |

`A`/`C`/`G`/`T` and `N` carry the real Phred quality character from `samtools`; `D`, `I` and `U` are placeholders and are written with quality `!` (Q0).

By default a UMI family is discarded outright as soon as any position resolves to `N`, `D`, `I` or `U`, so only families with a complete set of clean substitution calls reach the haplotype tables. Pass `--keep_nonconsensus` to retain those families with the codes above in place — useful for diagnosing whether a target is losing families to poor coverage (`U`), indel-prone alignment (`D`/`I`), or genuine ambiguity (`N`).

#### Why indels are discarded

The sites in `--snplist` were chosen on the assumption that **none of them should carry an indel**. That assumption comes from checking the candidate positions against the HGDP+1KG genomic data, where every one of them varies only by substitution — hence the name of the bundled list, [nonpolysnp_nonindel.csv](data/nonpolysnp_nonindel.csv). Given that, a `D` or `I` at one of these positions is treated as evidence that the UMI family is untrustworthy — typically an ONT indel error or a local misalignment — rather than as biology, so the whole minihaplotype is thrown out.

This is a population-level assumption, not a guarantee. A rare or previously unobserved sample may genuinely carry an indel at one of these sites, and the pipeline as it stands would silently discard exactly the families carrying it: the true variant never reaches `_haps.csv`, and the sample would instead be represented only by whichever families happen to look like clean substitutions. Nothing in the output flags that this has happened.

If you are working with a sample that might be outside the reference panel's diversity, re-run with `--keep_nonconsensus` and inspect `<barcode>_consensustable.csv`. A real indel should show up as `D` or `I` at a **consistent position across many independent UMI families**, whereas sequencing error scatters across positions and families. A target that loses a large and reproducible fraction of its families to the same `D`/`I` position is worth investigating rather than filtering away.

## Testing

Tests use [nf-test](https://www.nf-test.com/):

```bash
nf-test test tests/main.nf.test
```

The test runs the full pipeline against the bundled fixtures in [tests/input/](tests/input/) and [tests/data/](tests/data/).

## Repository layout

```
main.nf                   entry workflow (merge → trim → map → consensus → publish)
generate_consensus_onwards.nf  alternate entry point starting from an existing BAM (--bam)
workflows/                sub-workflows (GENERATE_CONSENSUS)
modules/                  process definitions (MAP_READS, GROUP_UMI, ...)
bin/                      Python helpers (group.py, call_consensus.py, process_consensus.py)
Docker/                   Dockerfile + conda environment.yml
data/                     example BED / SNP list inputs
tests/                    nf-test suite and fixtures
nextflow.config           default params, profiles, manifest
nextflow_schema.json      parameter schema
```

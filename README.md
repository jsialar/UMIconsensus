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
| Call consensus | `CALL_CONSENSUS` | Python | Derive a consensus sequence for each UMI family. |
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
| `--keep_nonconsensus` | `false` | Keep consensus sequences containing `D`/`N`. |
| `--noise_threshold` | `0.05` | Drop haplotypes below this fraction of the most abundant one. |
| `--minimap2_param` | `-ax map-ont -k 13 --MD` | minimap2 alignment parameters. |
| `--additional_outputs` | `false` | Also output consensus FASTQ files and amplicon BAMs. |
| `--debug` | `false` | Keep the Nextflow work directory (skip cleanup). |

See [nextflow_schema.json](nextflow_schema.json) for the full schema.

## Output

Results are published under `<output>/<barcode>/<module>/`:

```
<output>/<barcode>/
├── trim/           cutadapt JSON (and untrimmed reads if enabled)
├── aligned/        sorted BAM + index
├── group_umi/      UMI family size distributions (*_umisizes.tsv)
└── call_consensus/
    ├── *_consensustable.csv   consensus sequences per UMI family
    ├── *_quals.csv            per-base quality of dominant haplotypes
    └── *_haps.csv             dominant haplotype sequences
```

## Testing

Tests use [nf-test](https://www.nf-test.com/):

```bash
nf-test test tests/main.nf.test
```

The test runs the full pipeline against the bundled fixtures in [tests/input/](tests/input/) and [tests/data/](tests/data/).

## Repository layout

```
main.nf                   entry workflow (merge → trim → map → consensus → publish)
workflows/                sub-workflows (GENERATE_CONSENSUS)
modules/                  process definitions (MAP_READS, GROUP_UMI, ...)
bin/                      Python helpers (group.py, call_consensus.py, process_consensus.py)
Docker/                   Dockerfile + conda environment.yml
data/                     example BED / SNP list inputs
tests/                    nf-test suite and fixtures
nextflow.config           default params, profiles, manifest
nextflow_schema.json      parameter schema
```

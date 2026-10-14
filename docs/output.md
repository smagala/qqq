# cdcent/oamd-bio-qqq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Falco](#falco) - Raw read quality control
- [fastp](#fastp) - Short-read trimming and quality control
- [fastplong](#fastplong) - Long-read trimming and quality control
- [CAT_FASTQ](#cat_fastq) - Merge multiple runs of the same sample and platform
- [Dehosting](#dehosting) - Remove host reads using Deacon, Hostile, or Hostess
- [fadso](#fadso) - Optional read downsampling
- [MultiQC](#multiqc) - Aggregate quality control report
- [Pipeline information](#pipeline-information) - Workflow execution reports

### Falco

<details markdown="1">
<summary>Output files</summary>

- `falco/`
  - `*_fastqc_report.html`: FastQC-style report containing quality metrics, one per input FASTQ.
  - `*_fastqc_data.txt`: Quality metrics as tab-delimited text, one per input FASTQ.
  - `*_summary.txt`: Per-module pass/warn/fail summary, one per input FASTQ.

</details>

[Falco](https://falco.readthedocs.io/) is a drop-in C++ reimplementation of FastQC. It gives general quality metrics about your sequenced reads: quality score distribution, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. These metrics are computed on the **raw** reads, before dehosting.

The `_fastqc_data.txt` files are what MultiQC reads. MultiQC has no separate Falco module — its FastQC module parses these files, recognises the `##Falco` header, and labels the report section "Falco".

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*.fastp.fastq.gz`: Trimmed and filtered short reads.
  - `*.fastp.json`: Quality control and trimming statistics.
  - `*.fastp.html`: HTML quality control report.
  - `*.fastp.log`: Processing log.

</details>

[fastp](https://github.com/OpenGene/fastp) performs quality filtering and adapter trimming on short-read sequencing data, including single-end and paired-end reads.

Each sequencing run is processed separately before merging. This step is enabled by default and can be skipped using `--skip_fastp`.

### fastplong

<details markdown="1">
<summary>Output files</summary>

- `fastplong/`
  - `*.fastplong.fastq.gz`: Trimmed and filtered long reads.
  - `*.fastplong.json`: Quality control and trimming statistics.
  - `*.fastplong.html`: HTML quality control report.
  - `*.fastplong.log`: Processing log.

</details>

[fastplong](https://github.com/OpenGene/fastplong) performs quality filtering and trimming on long-read sequencing data, including Nanopore and PacBio reads.

Each sequencing run is processed separately before merging. This step is enabled by default and can be skipped using `--skip_fastp`.

### CAT_FASTQ

<details markdown="1">
<summary>Output files</summary>

- `cat/`
  - `*.merged.fastq.gz`: Merged FASTQ files from multiple sequencing runs.

</details>

CAT_FASTQ combines multiple sequencing runs with the same sample ID and sequencing platform.

Single-run samples bypass this step. Runs from different sequencing platforms are processed separately, even when they share the same sample ID.

For paired-end reads, R1 and R2 are merged separately.

Merging occurs after trimming when fastp or fastplong is enabled. When `--skip_fastp` is used, the original reads are merged instead.


### Dehosting

The pipeline supports three dehosting tools: Deacon, Hostile, and Hostess. The tool is selected using `--dehoster`.

Deacon is the default. Hostile and Hostess use Bowtie2 for short reads and minimap2 for long reads.

Dehosting can be disabled using `--dehoster none`.

#### Deacon

<details markdown="1">
<summary>Output files</summary>

- `deacon/`
  - `*.clean*.fq.gz`: Dehosted reads.
  - `*.deacon.json`: Dehosting statistics.

</details>

[Deacon](https://github.com/bede/deacon) removes potential host reads using minimizer-based sequence filtering.

#### Hostile

<details markdown="1">
<summary>Output files</summary>

- `hostile/`
  - `*.clean*.fastq.gz`: Dehosted reads.
  - `*.json`: Dehosting statistics.

</details>

[Hostile](https://github.com/bede/hostile) removes host reads using sequence alignment against a host reference.

#### Hostess

<details markdown="1">
<summary>Output files</summary>

- `hostess/`
  - `*.clean*.fastq.gz`: Dehosted reads.
  - `*.hostess.json`: Dehosting statistics.
  - `*.host.bam`: Host-aligned reads, when `--host_bam` is enabled.
  - `*.host.bam.bai`: BAM index, when host BAM output is enabled.
  - `*.host.stats.txt`: Alignment statistics, when generated.
  - `*.host.flagstat.txt`: Alignment summary, when generated.

</details>

Hostess performs alignment-based host removal and supports retaining host-aligned reads as BAM files when `--host_bam` is enabled.

Dehosting statistics are collected into the MultiQC report.

### fadso

<details markdown="1">
<summary>Output files</summary>

- `fadso/`
  - `*.downsampled.fastq.gz`: Downsampled reads.

</details>

[fadso](https://github.com/hkuwahara/sleaping) caps per-sample depth using s-leaping approximate reservoir sampling. This step only runs when `--max_reads` is set, and only for samples that exceed the cap — samples below it bypass the step entirely and appear unchanged in the pipeline's output reads. Paired-end mates are downsampled in lockstep so they stay synchronised.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

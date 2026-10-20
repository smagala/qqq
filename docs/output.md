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
- [BBDuk](#bbduk) - Optional ribosomal RNA depletion
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

### BBDuk

<details markdown="1">
<summary>Output files</summary>

- `bbduk/`
  - `*.depleted*.fastq.gz`: Reads that shared no k-mer with the rRNA reference. These continue through the pipeline.
  - `*.matched*.fastq.gz`: Reads that did. These are a deliverable, not a byproduct -- they are the input to rRNA characterisation, and they are the only way an over-aggressive reference shows up as something other than a quietly smaller library.
  - `*.bbduk.txt`: Per-reference-record hit counts, with the input and matched totals in the header. The complete file.
  - `*.bbduk.summary.txt`: The same header plus the top 20 records only. **This is what MultiQC reads**, not the complete file -- see below.
  - `*.bbduk.rpkm.txt`: Per-reference-record length, bases, coverage, read count and RPKM. Prefer this over `*.bbduk.txt` for anything quantitative: it carries the record length, without which read counts over-weight long references, since 23S is roughly twice 16S.
  - `*.bbduk.refstats.txt`: One row per reference *file*, so with a single `--rrna_ref` it is a one-line summary of reads and bases matched.

</details>

[BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) splits reads on whether they share a k-mer with a reference. This step runs only when `--rrna_ref` points at an rRNA FASTA, and sits after dehosting and before downsampling: before dehosting it would spend its pass on host reads that are about to be discarded, and after downsampling a `--max_reads` cap would already have been spent on ribosomes.

k-mer matching rather than alignment is the point. On a reference of tens of thousands of near-identical rRNA records a minimizer aligner's occurrence cap discards exactly the seeds that would identify the right record. Measured on a 19.4M-fragment wastewater RNA library against a 29,107-record RefSeq targeted-loci panel: an aligner found 24.5% of the library at identity >= 0.90, blastn found 68.4% at the identical threshold, and bbduk at `k=25` found 94.1% -- at 99.4% sensitivity on known rRNA fragments and 0.07% false positives on non-ribosomal bacterial 150-mers. `--rrna_k` is the only setting worth touching: `k=31` costs 2.4-3.1 points of sensitivity, `k=21` buys 0.2 points for 2.6x the false positives.

#### The reference

The reference is a versioned build-time asset, not something to assemble per run. The one AMDP publishes is at `s3://amdp-prd-analysis-static-assets/bbduk/silva138.2-ssu-lsu-nr99_rfam15.1-5s_20261020/`: SILVA 138.2 SSU+LSU Ref NR99 truncated exports plus Rfam 15.1 5S, U->T converted, 1,200,640 records. Its `README.md` is the authority on how it was built, what was deliberately left out, and what the per-record outputs may be read as. On an un-depleted total-RNA wastewater library it removes **95.75%** at 0.16% false positives.

Calibrate the expected removal rate to the library prep rather than to a single number: roughly 95% for un-depleted total RNA, single digits to ~15% for a bench-depleted RNA library, and ~0.5% for DNA, where the rRNA operon is a small fraction of a genome. A rate far outside the range for the stated prep is a sample-identity signal, not a tuning problem.

#### Reading the per-record output

Because SILVA deflines carry the full lineage, each row of `*.bbduk.rpkm.txt` reads as taxonomy. **Treat it as provenance, not as a composition profile.** bbduk credits each read to its best-matching record and breaks ties on reference file order, and on a reference where the right answer is often "one of these 400 near-identical records" the tie-break does most of the work. Reversing the order of all 1,200,640 records -- same reads, same reference content -- moves the composition by 0.25% at domain, 3.1% at phylum, 35.3% at genus and 71.7% at species; a different ambiguity policy moves phylum by 11.9%. Only the domain split is safe to quote. Per-record *coverage* stays meaningful regardless of which record earned the credit.

Depletion itself is unaffected by any of this, because it only asks whether any k-mer matched: reversing all 1.2M records changes the matched count by zero reads. Genus- and species-level claims need a classifier with LCA run on the `*.matched*.fastq.gz` stream, which is why those reads are published.

One parsing caveat: splitting a SILVA lineage on `;` and taking field *n* is right for Bacteria and Archaea, which are six ranks then an organism name, but wrong for Eukaryota, whose lineages are much deeper. Use SILVA's own rank map (`taxonomy/tax_slv_{ssu,lsu}_138.2.txt.gz`) if you need ranks.

#### Why MultiQC gets the truncated file

MultiQC writes every row of the file it parses into its own outputs. `*.bbduk.txt` has one row per matched reference record, and against a reference large enough to be worth using that is a lot of rows -- 93,601 on the validation library, which cost 13 MB of `bbmap.txt` and 23 MB of `multiqc_data.json` for a **single sample**. A hundred-sample run would be a gigabyte of JSON. `*.bbduk.summary.txt` parses identically, because MultiQC's `bbmap/stats` pattern matches on the file header rather than its name, and yields the same percent-filtered figure from 38 KB.

bbduk's stderr log is deliberately not wired into MultiQC at all: its `bbduk` pattern looks for `Executing jgi.BBDuk`, which BBTools renamed to `bbduk.BBDukS` in 39.76, and which the JVM puts on line three where the pattern allows two. It matches no version this pipeline runs.

A counts footnote: `#Matched` in the stats file runs 0.27-0.28% below the `Contaminants` figure bbduk reports to its log. The split is pair-aware -- one matching mate removes both -- while attribution is per-mate evidence. Neither is wrong; they answer different questions.

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

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
- [rasusa](#rasusa) - Optional read downsampling
- [MultiQC](#multiqc) - Aggregate quality control report
- [Run report](#run-report) - What processed the run, what each step removed, and where every read went
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
  - `*.paired.fail.fastq.gz`, `*_R1.fail.fastq.gz`, `*_R2.fail.fastq.gz`: Reads the filter removed (paired-end), when `--save_failed_reads` is on (the default). `*.paired.fail.fastq.gz` holds each failing read, its name tagged with the reason (`failed_quality_filter`, `failed_too_short`, ...); `*_R1/_R2.fail.fastq.gz` hold their mates, which passed but were dropped with the pair. **Not complete:** fastp writes no file at all for a pair in which *both* mates fail, while still counting them as removed -- on the test sample 32 of 373 removed pairs. The run report marks these entries `"complete": false`.
  - `*.fail.fastq.gz`: Reads the filter removed (single-end), when `--save_failed_reads` is on. Complete.

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
  - `*.fail.fastq.gz`: Reads the filter removed, names tagged with the reason, when `--save_failed_reads` is on (the default). Complete.

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

The reference is a versioned build-time asset, not something to assemble per run. The one AMDP publishes, and the oamd-bio-nf-configs default, is at `s3://amdp-prd-analysis-static-assets/bbduk/silva138.2-ssu-lsu-nr99_rfam15.1-5s_mask20261023/`: SILVA 138.2 SSU+LSU Ref NR99 truncated exports plus Rfam 15.1 5S, U->T converted, 1,200,640 records, with 4,299 bp of audited non-ribosomal blocks in 10 records masked to `N` (`PANEL_CARD.md` beside it). The unmasked 2026-10-20 build stays at `.../silva138.2-ssu-lsu-nr99_rfam15.1-5s_20261020/`, and its `README.md` is the authority on how the base was built, what was deliberately left out, and what the per-record outputs may be read as. On an un-depleted total-RNA wastewater library the unmasked build removes **95.75%** at 0.16% false positives.

Calibrate the expected removal rate to the library prep rather than to a single number: roughly 95% for un-depleted total RNA, single digits to ~15% for a bench-depleted RNA library, and ~0.5% for DNA, where the rRNA operon is a small fraction of a genome. A rate far outside the range for the stated prep is a sample-identity signal, not a tuning problem.

#### Reading the per-record output

Because SILVA deflines carry the full lineage, each row of `*.bbduk.rpkm.txt` reads as taxonomy. **Treat it as provenance, not as a composition profile.** bbduk credits each read to its best-matching record and breaks ties on reference file order, and on a reference where the right answer is often "one of these 400 near-identical records" the tie-break does most of the work. Reversing the order of all 1,200,640 records -- same reads, same reference content -- moves the composition by 0.25% at domain, 3.1% at phylum, 35.3% at genus and 71.7% at species; a different ambiguity policy moves phylum by 11.9%. Only the domain split is safe to quote. Per-record *coverage* stays meaningful regardless of which record earned the credit.

Depletion itself is unaffected by any of this, because it only asks whether any k-mer matched: reversing all 1.2M records changes the matched count by zero reads. Genus- and species-level claims need a classifier with LCA run on the `*.matched*.fastq.gz` stream, which is why those reads are published.

One parsing caveat: splitting a SILVA lineage on `;` and taking field *n* is right for Bacteria and Archaea, which are six ranks then an organism name, but wrong for Eukaryota, whose lineages are much deeper. Use SILVA's own rank map (`taxonomy/tax_slv_{ssu,lsu}_138.2.txt.gz`) if you need ranks.

#### Why MultiQC gets the truncated file

MultiQC writes every row of the file it parses into its own outputs. `*.bbduk.txt` has one row per matched reference record, and against a reference large enough to be worth using that is a lot of rows -- 93,601 on the validation library, which cost 13 MB of `bbmap.txt` and 23 MB of `multiqc_data.json` for a **single sample**. A hundred-sample run would be a gigabyte of JSON. `*.bbduk.summary.txt` parses identically, because MultiQC's `bbmap/stats` pattern matches on the file header rather than its name, and yields the same percent-filtered figure from 38 KB.

bbduk's stderr log is deliberately not wired into MultiQC at all: its `bbduk` pattern looks for `Executing jgi.BBDuk`, which BBTools renamed to `bbduk.BBDukS` in 39.76, and which the JVM puts on line three where the pattern allows two. It matches no version this pipeline runs.

A counts footnote: `#Matched` in the stats file runs 0.27-0.28% below the `Contaminants` figure bbduk reports to its log. The split is pair-aware -- one matching mate removes both -- while attribution is per-mate evidence. Neither is wrong; they answer different questions.

### rasusa

<details markdown="1">
<summary>Output files</summary>

- `rasusa/`
  - `*.downsampled.fastq.gz`: Downsampled reads.

</details>

[rasusa](https://github.com/mbhall88/rasusa) caps per-sample depth. `rasusa reads --num` measures the input before sampling, so the cap is an exact count rather than an approximation, and paired mates are kept or dropped as templates so they stay synchronised. This step only runs when `--max_reads` is set, and only for samples that exceed the cap — samples below it bypass the step entirely and appear unchanged in the pipeline's output reads.

The bypass is an optimisation, not a correctness requirement: an input holding fewer reads than the cap is written out whole, with exit 0. Rerunning with `--dehoster none` routes every sample through the module instead, because there is no dehoster read count to gate on, and the samples under the cap come out of it unchanged.

Reproducing a downsampled run needs the rasusa version as well as `--downsample_seed`. Which reads a given seed selects changed between rasusa major versions, so the seed alone does not identify the selection; the version is recorded in the software versions section of the MultiQC report.

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

### Run report

<details markdown="1">
<summary>Output files</summary>

- `run_report/`
  - `qqq_run_report.html`: Self-contained page for people: headline tiles, a sortable and filterable per-sample table with a bar showing where each sample's input went, the tools, versions and arguments, the databases, and the published location of every sample's reads.
  - `qqq_run_report.json`: The same information for machines, with the schema below.

</details>

Written once, after MultiQC, by `QQQ_RUN_REPORT` (`bin/qqq_run_report.py`). It answers three questions that otherwise mean reading a dozen directories: exactly what processed the run, what each step removed from each sample, and where the reads -- kept and removed -- are.

#### JSON schema (`schema_version` 1.0.0)

The schema is versioned by `schema_version`, [semver](https://semver.org/): a new field is a minor bump, a removed or re-typed field a major one. Consumers should check the major version.

**Units.** Read counts are sequence records; paired-end mates are counted separately. That is how every tool in the chain counts (fastp `total_reads`, deacon `seqs_in`, hostile `reads_in`, bbduk `#Total`), so each step's input equals the previous step's output without halving anything. Counts that could not be measured -- a stub run, a missing report -- are `null`, never `0`.

**Paths.** Every published file is an object `{"path": ..., "relative_path": ...}`. `path` is under `params.outdir` exactly as given when that is a URI (`s3://...`), or made absolute when it is local; `relative_path` is relative to the outdir, for consumers that have moved or synced it.

Top level:

| Field | Type | Meaning |
| ----- | ---- | ------- |
| `schema_version` | string | Version of this schema. |
| `generated_at` | string | ISO 8601 UTC time the report was written. |
| `generator` | object | `script` and the `python` version that wrote it. |
| `units` | string | The counting convention above, in words. |
| `step_order` | array | `["read_filter", "host_removal", "rrna_depletion", "downsampling"]`: the removal steps in the order reads pass through them. |
| `run` | object | `pipeline`, `pipeline_version`, `commit` (from `workflow.commitId`, or `git rev-parse HEAD` in the project directory when the run was launched from a checkout), `commit_source`, `commit_dirty` (tracked files differed from that commit; `null` when unknown), `revision`, `repository`, `run_name`, `session_id`, `started`, `nextflow_version`, `profile`, `container_engine`, `outdir`, `launch_dir`, `work_dir`, `project_dir`, `command_line`, and `params`: the parameters that decide what the run did. |
| `tools` | array | One entry per tool: `name`, `version` (`null` when the module's version command returned nothing, with `version_note`), `versions` (all distinct), `processes`, `arguments`, `arguments_source`. |
| `databases` | array | One entry per database or reference: `role` (`host_index`, `host_reference_fasta`, `rrna_reference`), `tool`, `path` (full), `name` (file or directory name), `is_directory`, `size_bytes`, `md5`, `sha256`, `checksum_source`, `checksum_note`, `parameters`. |
| `summary` | object | Run totals: `samples`, `samples_zero_reads`, `samples_zero_reads_ids`, `samples_with_warnings`, `input_reads`, `input_bases`, `output_reads`, `retained_fraction`, `removed_reads` and `removed_fraction` keyed by step, `read_filter_removed_by_reason`, `adapter_trimmed_reads`, `adapter_trimmed_fraction`, `adapter_trimmed_bases`. |
| `outputs` | object | Run-level files: `run_report_json`, `run_report_html`, `multiqc_report`, `software_versions`. |
| `samples` | array | One entry per sample, below. |

**Arguments.** `arguments_source` says where an entry came from. `tool_report` means the tool recorded it at runtime: the fastp/fastplong `command` (per-sample file names stripped, one entry per distinct set with a `samples` count), deacon's `k`, `w`, `abs_threshold`, `rel_threshold` and switches, hostile/hostess's aligner, index and options. `pipeline_config` means the pipeline's own `conf/modules.config` or a module input -- bbduk's `k=<rrna_k> interleaved=f`, rasusa's `--num`/`--seed`, hostess's `--bam` -- so a launch-time `-c` override of `ext.args` would not show there.

**Checksums.** Nothing large is read over the network to produce them: the deacon index and the rRNA reference are GB on S3. They come from a sidecar beside the file (`<file>.md5`, `.md5sum`, `.sha256`, `.sha256sum`, or a `md5sums.txt`/`MD5SUMS`/`sha256sums.txt`/`SHA256SUMS` in the same directory naming the file), or are computed for local files of 64 MiB or less. Otherwise they are `null` and `checksum_note` says why. For a host database the tool's own runtime view (`parameters.reported_by_tool`) is attached: for deacon, the `k`, `w` and thresholds from its summary.

Each `samples[]` entry:

| Field | Type | Meaning |
| ----- | ---- | ------- |
| `id` | string | `<sample>_<platform>`, the identity the sample's files are named by. |
| `single_end`, `platform`, `read_type` | | As routed (`read_type` is `short` or `long`). |
| `runs` | array | `{id, input_files}` per samplesheet run, with the input FASTQ locations. |
| `input_reads`, `input_bases` | integer | What entered the first step that ran. |
| `output_reads` | integer | What left the last step: the final read count. |
| `retained_fraction` | number | `output_reads / input_reads`. |
| `zero_reads` | boolean | No reads left. |
| `steps` | object | Keyed by every name in `step_order`, always all four. |
| `outputs` | object | `clean_reads` and `removed_reads`, below. |
| `warnings` | array | Strings: a step whose input does not match the previous step's output, a report that could not be read, zero reads left. |

Each `steps.<step>` has `status` (`ran`, `skipped` -- off for the run, with `reason` -- or, for downsampling only, `bypassed`: the sample was under `max_reads` and never reached rasusa), `tool`, `reason`, `input_reads`, `removed_reads` and `output_reads`. Steps add:

- `read_filter`: `input_bases`, `output_bases`, `removed_by_reason` (fastp's `filtering_result`: `low_quality`, `too_many_N`, `too_short`, ...), `adapter_trimmed_reads`, `adapter_trimmed_bases`, `adapter_trimming_effective` (whether adapter removal cut anything at all), and `runs`, the per-run figures and full fastp command. Summed over runs for a merged sample.
- `host_removal`: `input_bases`, `output_bases`, `removed_bases` (deacon only; hostile and hostess do not report bases), `report`.
- `rrna_depletion`: from the bbduk stats header, `#Total` in and `#Matched` removed; `report`. For paired-end samples it also carries `exact: false` and a `note`: `#Matched` counts reads whose own k-mers matched, but bbduk removes the whole pair, so the records actually removed run slightly higher -- 216 against 210 on the test sample, 0.27-0.28% on the validation library (see the counts footnote under [BBDuk](#bbduk)). The pair-aware figure is printed only to bbduk's stderr, which the module does not keep.
- `downsampling`: `cap_reads` (max_reads in records), and `derived: true`. rasusa writes no report; its counts follow from `rasusa reads --num` keeping exactly the cap, or everything when the input is under it, which the module's nf-test pins.

`outputs.clean_reads` is `{step, files, reason}`: the step that last wrote the sample's reads and the published files. When no step rewrote them (a single-run sample with every step off) `files` is `null`, `reason` says so, and `source_files` gives the samplesheet inputs.

`outputs.removed_reads.<step>` is `{files, reason}` for every step in `step_order`: either the published files and `reason: null`, or `files: null` and a `reason`:

| Step | Saved | Otherwise |
| ---- | ----- | --------- |
| `read_filter` | fastp/fastplong `*.fail.fastq.gz`, all runs, with `--save_failed_reads` (default). Paired-end fastp entries carry `complete: false` and a `note`: see [fastp](#fastp). | `save_failed_reads=false; ...`, or `step skipped (skip_fastp=true)` |
| `host_removal` | hostess `*.host.bam` and `.bai`, with `--host_bam` | `host reads not retained (deacon/hostile keep none)`, `host_bam=false; ...`, or `step skipped (dehoster=none)` |
| `rrna_depletion` | bbduk `*.matched*.fastq.gz` | `step skipped (rrna_ref not set)` |
| `downsampling` | never written | rasusa writes only what it keeps; the complete pre-downsampling reads are the previous step's published output |

An abbreviated example:

```json
{
  "schema_version": "1.0.0",
  "run": { "pipeline_version": "1.0.0dev", "commit": "6ae8454d...", "run_name": "pedantic_feynman", "outdir": "s3://bucket/results/run1" },
  "tools": [
    { "name": "deacon", "version": "0.17.0", "processes": ["DEACON_FILTER"], "arguments_source": "tool_report",
      "arguments": [{ "index": "host.k31w15.idx", "settings": { "k": 31, "w": 15, "abs_threshold": 2, "rel_threshold": 0.01, "ordered": true }, "samples": 48 }] }
  ],
  "databases": [
    { "role": "rrna_reference", "tool": "bbduk", "path": "s3://.../silva138.2-ssu-lsu-nr99_rfam15.1-5s_mask20261023.fa.gz",
      "name": "silva138.2-ssu-lsu-nr99_rfam15.1-5s_mask20261023.fa.gz", "md5": null, "sha256": null, "parameters": { "k": 25 } }
  ],
  "summary": { "samples": 48, "input_reads": 951203110, "output_reads": 40211876, "retained_fraction": 0.0423, "samples_zero_reads": 1 },
  "samples": [
    {
      "id": "S1_illumina", "single_end": false, "input_reads": 19860772, "output_reads": 63652, "zero_reads": false,
      "steps": {
        "read_filter":    { "status": "ran", "tool": "fastp",  "input_reads": 19860772, "removed_reads": 957480,  "output_reads": 18903292, "adapter_trimmed_reads": 11859090 },
        "host_removal":   { "status": "ran", "tool": "deacon", "input_reads": 18903292, "removed_reads": 9869852, "output_reads": 9033440 },
        "rrna_depletion": { "status": "ran", "tool": "bbduk",  "input_reads": 9033440,  "removed_reads": 8969788, "output_reads": 63652 },
        "downsampling":   { "status": "skipped", "reason": "max_reads not set", "input_reads": null, "removed_reads": null, "output_reads": null }
      },
      "outputs": {
        "clean_reads": { "step": "rrna_depletion", "files": [{ "path": "s3://bucket/results/run1/bbduk/S1_illumina.depleted_1.fastq.gz", "relative_path": "bbduk/S1_illumina.depleted_1.fastq.gz" }, "..."] },
        "removed_reads": {
          "read_filter":    { "files": [{ "path": "s3://bucket/results/run1/fastp/S1_illumina_run1.paired.fail.fastq.gz", "relative_path": "..." }, "..."], "reason": null, "complete": false },
          "host_removal":   { "files": null, "reason": "host reads not retained (deacon/hostile keep none)" },
          "rrna_depletion": { "files": [{ "path": "s3://bucket/results/run1/bbduk/S1_illumina.matched_1.fastq.gz", "relative_path": "..." }, "..."], "reason": null },
          "downsampling":   { "files": null, "reason": "step skipped (max_reads not set)" }
        }
      },
      "warnings": []
    }
  ]
}
```

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

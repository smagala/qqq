# cdcent/oamd-bio-qqq: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

QQQ processes short- and long-read sequencing data, including Illumina, Nanopore, and PacBio. The pipeline uses Falco for raw read QC, fastp or fastplong for optional trimming, and combines multiple runs from the same sample and platform. It also supports host read removal and optional downsampling. MultiQC is used to summarize the results.

### Running on AMDP

Pipeline configuration (container tags, the dehosting index locations, resource labels)
lives in [oamd-bio-nf-configs](https://github.com/cdcent/oamd-bio-nf-configs) rather
than in this repository, and is loaded at runtime:

```bash
export NF_PIPELINE_REPO=oamd-bio-qqq          # selects conf/pipeline/oamd-bio-qqq/
export NF_CONFIGS_REPO=/abs/path/to/oamd-bio-nf-configs
export AWS_PROFILE=<profile for the target tier>

nextflow -c "$NF_CONFIGS_REPO/conf/pipeline/custom.config" \
  run cdcent/oamd-bio-qqq \
  -profile prdanalysis \
  --input samplesheet.csv \
  --outdir s3://<bucket>/<prefix>
```

Three things about those variables are easy to get wrong, and none of them fail
loudly on their own:

- **`NF_PIPELINE_REPO` has no safe default.** Unset, it falls back to
  `oamd-bio-ph-core`, which silently selects that pipeline's configuration and
  sends `workDir` to `s3://amdp-<tier>-bio-development-work/oamd-bio-ph-core/`.
- **`NF_CONFIGS_REPO` must be an absolute path.** It is interpolated into
  `includeConfig`, which resolves relative paths against the config file's own
  directory rather than your working directory, producing a doubled path.
- **The checkout it points at must contain `conf/pipeline/oamd-bio-qqq/`.** If it
  does not, the configs are skipped without error and the run reaches the dehoster
  with no reference index.

`run-qqq-prda.sh` in the parent directory of a standard checkout sets all three
correctly for the prd-analysis tier.

### Read trimming and merging

By default, QQQ performs read trimming and quality filtering after raw read QC.

- **fastp** is used for short-read sequencing data.
- **fastplong** is used for long-read sequencing data.
- `--skip_fastp` disables trimming for both short and long reads.

Each sequencing run is processed separately before merging.

When multiple runs share the same sample ID and sequencing platform, CAT_FASTQ combines them into a single set of FASTQ files.

- Single-run samples bypass merging.
- Runs from different sequencing platforms are processed separately.
- Paired-end reads are merged separately for R1 and R2.

When trimming is enabled, merging occurs after fastp or fastplong. When `--skip_fastp` is enabled, the original FASTQ files are merged instead.

Trimming results are saved in the `fastp/` and `fastplong/` output directories. Merged reads are saved in `cat/`.

### Choosing a dehoster

`--dehoster` selects the host removal implementation, or turns the step off.

| `--dehoster` | Method | Host BAM available? | Notes |
| ------------ | ------ | ------------------- | ----- |
| `deacon` | Minimizer-based filtering | No | Default; fast, approximate host filtering |
| `hostess` | Bowtie2 (short reads) / minimap2 (long reads) | Yes, with `--host_bam` | Alignment-based host removal |
| `hostile` | Bowtie2 (short reads) / minimap2 (long reads) | No | Alignment-based host removal |
| `none` | No dehosting | No | Skips host removal; QC, trimming, and merging remain available |

`deacon` does not align. It counts how many of each read's minimizers appear in a
prebuilt index and drops the read when enough of them do, which is why it is so
much faster than the other two — and why it is approximate where they are exact.

Its output is BGZF when the image carries `bgzip`, as the AMDP image does, and
plain gzip otherwise; both are valid gzip. `--ordered` is passed by default so a
rerun reproduces the same bytes, which deacon does not otherwise guarantee.

Two consequences are worth knowing before switching:

- **`--host_bam` has no effect under `deacon`.** There is no alignment to keep.
  Use `--dehoster hostess` if you need the removed reads as a BAM.
- **The reference matters, and it is yours to choose.** The pipeline ships no
  default index: `--deacon_db` and `--hostile_db` are unset unless a deployment
  config supplies them. Point every dehoster at builds of the *same* reference
  unless you have a reason not to, or `% Host` will not be comparable between
  them. AMDP deployments do exactly that, via `oamd-bio-nf-configs`.
- **`deacon` will not agree with the aligners exactly, even on one reference.**
  It is approximate by design. Measured across the 49-sample talos set with all
  three reading the same masked reference, deacon kept about 0.5% fewer reads
  than `hostile`/`hostess`, and 0.36% fewer on the worst single sample. Reading
  `panhuman-1` instead — a *different* reference, masked less aggressively —
  roughly doubled that gap, to 0.7–0.95% typical and 1.54% worst-case.
  `hostile` and `hostess` agree with each other read for read.

### Dehosting reference

For `deacon`, `--deacon_db` must point to a prebuilt minimizer index file (`*.idx`).

For `hostile` and `hostess`, `--hostile_db` must point to the directory containing the required reference indexes.

- Short reads use Bowtie2 indexes (`<index>.*.bt2`).
- Long reads use a minimap2 index (`<index>.mmi`).
- `--hostile_index` specifies the index basename.

When both short-read and long-read data are processed, the reference directory must contain the indexes required for both alignment methods.

Alternatively, `--host_fasta` can be used to build the required indexes from a reference FASTA during pipeline execution. When provided, it takes precedence over `--deacon_db` and `--hostile_db`.

### Downsampling

`--max_reads` caps how many reads (single-end) or read pairs (paired-end) each
sample keeps after dehosting; `--downsample_seed` (default `42`) makes the
sampling reproducible. Only samples above the cap are downsampled, and the
decision reuses the read counts the dehoster already reports. Leave it unset to
disable the step entirely.


## Samplesheet input

The pipeline accepts a CSV samplesheet containing information about the sequencing data.

The supported columns are `sample`, `fastq_1`, `fastq_2`, and `platform`.

The `platform` column is optional. If it is not provided, the pipeline defaults to `illumina`.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet columns

| Column | Description |
| ------ | ----------- |
| `sample` | Sample ID. Multiple sequencing runs from the same sample should use the same ID. |
| `fastq_1` | Required FASTQ file for read 1 or a single-end/long-read run. |
| `fastq_2` | FASTQ file for read 2 in paired-end sequencing. Leave empty for single-end and long-read data. |
| `platform` | Optional sequencing platform. Defaults to `illumina` if not provided. |

FASTQ input files must be gzipped and use the `.fastq.gz` or `.fq.gz` extension.

### Supported sequencing platforms

The pipeline supports both short-read and long-read platforms.

**Short-read platforms:** `illumina`, `ls454`, `abi_solid`, `bgiseq`, `complete_genomics`, `dnbseq`, `element`, `genapsys`, `genemind`, `helicos`, `ion_torrent`, `tapestri`, `ultima`, `vela_diagnostics`, and `capillary`.

**Long-read platforms:** `nanopore`, `ont`, `oxford_nanopore`, `pacbio`, and `pacbio_smrt`.

Short reads are processed with fastp, while long reads are processed with fastplong.

Each long-read run must contain a single FASTQ file in `fastq_1`, with `fastq_2` left empty.

### Multiple runs of the same sample

When a sample has multiple sequencing runs, use the same sample ID and platform for each run.

For example:

```csv
sample,fastq_1,fastq_2,platform
SAMPLE1,SAMPLE1_L001_R1.fastq.gz,SAMPLE1_L001_R2.fastq.gz,illumina
SAMPLE1,SAMPLE1_L002_R1.fastq.gz,SAMPLE1_L002_R2.fastq.gz,illumina
SAMPLE1,SAMPLE1_nanopore_run1.fastq.gz,,nanopore
SAMPLE1,SAMPLE1_nanopore_run2.fastq.gz,,nanopore
SAMPLE2,SAMPLE2_pacbio.fastq.gz,,pacbio
```

The pipeline processes each run separately for raw read QC and optional trimming.

After trimming, runs with the same sample ID and sequencing platform are combined using CAT_FASTQ.

Runs from different platforms are not merged together, even when they share the same sample ID.

Single-run samples bypass CAT_FASTQ.

When `--skip_fastp` is enabled, reads are merged without trimming.

### Example samplesheet

An [example samplesheet](../assets/samplesheet.csv) is provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run cdcent/oamd-bio-qqq --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run cdcent/oamd-bio-qqq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull cdcent/oamd-bio-qqq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [cdcent/oamd-bio-qqq releases page](https://github.com/cdcent/oamd-bio-qqq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

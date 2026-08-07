# cdcent/oamd-bio-qqq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of cdcent/oamd-bio-qqq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Optional read downsampling after dehosting with `fadso`, controlled by `--max_reads` and `--downsample_seed`. Only samples above the cap are routed to the module; the decision reuses the read counts Hostile already reports, so nothing makes an extra pass over the FASTQs. Disabled by default. ([AMD-31664](https://cdc-oamd.atlassian.net/browse/AMD-31664))
- `--hostile_fasta`, which builds the Hostile index at runtime instead of using a prebuilt one. Used by the test profiles so they do not depend on hosted binaries. ([AMD-31686](https://cdc-oamd.atlassian.net/browse/AMD-31686))
- Falco and Hostile citations in `CITATIONS.md`, replacing FastQC.

### `Fixed`

- **Falco results were silently missing from every MultiQC report.** MultiQC has no Falco module: its FastQC module parses Falco's output and renames the section on seeing the `##Falco` header. The `fastqc/data` search pattern was `fastqc_data.txt` (an exact filename) until MultiQC 1.32, and only then became `*fastqc_data.txt`. Falco emits `<name>_fastqc_data.txt`, so the pinned MultiQC 1.31 matched nothing. ([AMD-31663](https://cdc-oamd.atlassian.net/browse/AMD-31663))
- Falco paired-end output no longer overwrites itself. Falco writes to fixed output names, so a single invocation over both mates replaced R1's results with R2's. Each mate is now processed separately with names derived from the input file.
- `-profile test` and `-profile test_full` failed immediately with `Argument of file() function cannot be null`. Both were untouched template boilerplate pointing at another pipeline's test data, and neither set the mandatory `hostile_db`. ([AMD-31686](https://cdc-oamd.atlassian.net/browse/AMD-31686))
- The pipeline identity was `cdcent/qqq` throughout, including the manifest, `modules.json`, MultiQC section IDs and report links. It is now `cdcent/oamd-bio-qqq`.
- `nf-test` CI could not run at all: the composite action used a GitHub Action blocked by the enterprise allowed-actions policy, and that is enforced at workflow parse time, so every job failed within seconds. ([AMD-31693](https://cdc-oamd.atlassian.net/browse/AMD-31693))
- The sendmail template read a pipeline logo unconditionally. No logo ships with this pipeline, so `--email` raised `FileNotFoundException`. The image part is now emitted only when the file exists.
- Removed `ph_core_config_base` URL loading. `oamd-bio-nf-configs` is private, so those raw URLs returned 404 and broke every `nf-core` tooling command run against this repository. Pipeline configs are loaded from a local checkout instead.
- Removed unused iGenomes boilerplate (`params.fasta`, `getGenomeAttribute`) and the `hook_url` notification plumbing, whose upstream implementation no longer exists.

### `Dependencies`

| Tool           | Old version | New version |
| -------------- | ----------- | ----------- |
| MultiQC        | 1.31        | 1.35        |
| samtools       | 1.22.1_v9   | 1.22.1_v11  |
| falco          |             | 1.2.5       |
| hostile        |             | 2.0.2       |
| fadso          |             | 20230810    |
| bowtie2        |             | 2.5.4       |

- `falco` now comes from the public [amd-ph-core/modules](https://github.com/amd-ph-core/modules) repository rather than the private ph-core one.
- `hostile/clean`, `samtools/bam2fq` and `multiqc` updated to versions that emit tool versions through the `versions` topic channel. Every module in the pipeline now uses the topic; MultiQC deliberately does not, because its input depends on that topic resolving and emitting to it would deadlock the pipeline.
- `utils_nextflow_pipeline`, `utils_nfcore_pipeline` and `utils_nfschema_plugin` updated.

### `Deprecated`

- `--hook_url` and the associated Slack/Teams notification templates. The upstream `imNotification` helper was removed from nf-core and the feature no longer exists.

# cdcent/oamd-bio-qqq

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/cdcent/oamd-bio-qqq)
[![GitHub Actions CI Status](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/linting.yml/badge.svg)](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/cdcent/oamd-bio-qqq)

## Introduction

**cdcent/oamd-bio-qqq** (**Q**uick **Q**uality & **Q**uarantine) is a Nextflow pipeline for dehosting paired-end Illumina reads. It produces original read quality metrics, removes host reads, and generates an aggregated report. The pipeline is designed to be quick and easy to run across many samples that may contain host contamination.

Host removal is done by [Deacon](https://github.com/bede/deacon) by default, which matches read minimizers
against a prebuilt index rather than aligning, making it substantially faster than an aligner at the cost of
being approximate. `--dehoster hostess` and `--dehoster hostile` select bowtie2-based alternatives, and are
the ones to reach for when exactness matters or when the removed reads are wanted as a BAM; `--dehoster none`
turns the step off. See [usage](docs/usage.md#choosing-a-dehoster) for the trade-offs, including the fact that
each dehoster uses a differently masked reference and so will not remove exactly the same reads.

The pipeline performs the following steps:

1. Read QC ([`Falco`](https://falco.readthedocs.io/)) — quality metrics on raw reads before dehosting
2. Host read removal ([`Deacon`](https://github.com/bede/deacon), [`Hostile`](https://github.com/bede/hostile)) — decontaminate reads against a reference index
3. Read downsampling ([`fadso`](https://github.com/hkuwahara/sleaping)) — _optional_, cap the depth of the dehosted reads
4. Aggregate reports ([`MultiQC`](http://multiqc.info/)) — combined QC and dehosting summary

Downsampling is off by default. Setting `--max_reads` caps how many reads (single-end) or read pairs (paired-end) each sample keeps: samples above the cap are downsampled to it, and samples below it pass through untouched. The decision uses the read counts the dehoster already reports, so no extra pass over the FASTQs is needed. `--downsample_seed` (default `42`) keeps the sampling reproducible across runs.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with your input data:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
SAMPLE1,SAMPLE1_R1.fastq.gz,SAMPLE1_R2.fastq.gz
```

Each row represents a pair of paired-end FASTQ files.

Run the pipeline:

```bash
nextflow run cdcent/oamd-bio-qqq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --deacon_db /path/to/panhuman-1.k31w15.idx \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

cdcent/oamd-bio-qqq was written by James Smagala (OAMD Bioinformatics).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

This pipeline has not been released yet, so it has no DOI to cite. A Zenodo DOI will
be added here on first release.

References for the tools used by the pipeline are listed in [`CITATIONS.md`](CITATIONS.md).
The MultiQC report also carries a methods description naming the tools actually used
for that run.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

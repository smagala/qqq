# cdcent/oamd-bio-qqq

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/cdcent/oamd-bio-qqq)
[![GitHub Actions CI Status](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/linting.yml/badge.svg)](https://github.com/cdcent/oamd-bio-qqq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/cdcent/oamd-bio-qqq)

## Introduction

**cdcent/oamd-bio-qqq** (**Q**uick **Q**uality & **Q**uarantine) is a Nextflow pipeline for dehosting paired-end Illumina reads. It produces original read quality metrics, removes host reads, and generates an aggregated report. The pipeline is designed to be quick and easy to run across many samples that may contain host contamination.

By default, qqq uses a masked human reference genome that excludes any region potentially overlapping with pathogen sequences, ensuring pathogen reads are never inadvertently discarded during dehosting.

The pipeline performs the following steps:

1. Read QC ([`Falco`](https://falco.readthedocs.io/)) — quality metrics on raw reads before dehosting
2. Host read removal ([`Hostile`](https://github.com/bede/hostile)) — decontaminate reads against a reference index
3. Read downsampling ([`fadso`](https://github.com/hkuwahara/sleaping)) — _optional_, cap the depth of the dehosted reads
4. Aggregate reports ([`MultiQC`](http://multiqc.info/)) — combined QC and dehosting summary

Downsampling is off by default. Setting `--max_reads` caps how many reads (single-end) or read pairs (paired-end) each sample keeps: samples above the cap are downsampled to it, and samples below it pass through untouched. The decision uses the read counts Hostile already reports, so no extra pass over the FASTQs is needed. `--downsample_seed` (default `42`) keeps the sampling reproducible across runs.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

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
   --hostile_db /path/to/hostile/index \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

cdcent/oamd-bio-qqq was written by James Smagala (OAMD Bioinformatics).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use cdcent/oamd-bio-qqq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

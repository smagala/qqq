/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FALCO                  } from '../modules/amd-ph-core/falco/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTPLONG              } from '../modules/nf-core/fastplong/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FADSO                  } from '../modules/amd-ph-core/fadso/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { DEACON_INDEX           } from '../modules/amd-ph-core/deacon/index/main'
include { DEACON_FILTER as DEACON_FILTER_SMALL } from '../modules/amd-ph-core/deacon/filter/main'
include { DEACON_FILTER as DEACON_FILTER_LARGE } from '../modules/amd-ph-core/deacon/filter/main'
include { HOSTESS as HOSTESS_SMALL_SHORT } from '../modules/amd-ph-core/hostess/main'
include { HOSTESS as HOSTESS_LARGE_SHORT } from '../modules/amd-ph-core/hostess/main'
include { HOSTESS as HOSTESS_SMALL_LONG } from '../modules/amd-ph-core/hostess/main'
include { HOSTESS as HOSTESS_LARGE_LONG } from '../modules/amd-ph-core/hostess/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_SMALL_SHORT } from '../modules/nf-core/hostile/clean/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_LARGE_SHORT } from '../modules/nf-core/hostile/clean/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_SMALL_LONG } from '../modules/nf-core/hostile/clean/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_LARGE_LONG } from '../modules/nf-core/hostile/clean/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qqq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QQQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_multiqc_files = channel.empty()

    ch_reads = ch_samplesheet

    // Identify and branch reads by platform
    // (short-read vs long-read)
    //
    short_platforms = [
        "illumina",
        "ls454",
        "abi_solid",
        "bgiseq",
        "complete_genomics",
        "dnbseq",
        "element",
        "genapsys",
        "genemind",
        "helicos",
        "ion_torrent",
        "tapestri",
        "ultima",
        "vela_diagnostics",
        "capillary"
    ]

    long_platforms = [
        "nanopore",
        "ont",
        "pacbio",
        "pacbio_smrt",
        "oxford_nanopore"
    ]

    ch_reads
        .map { meta, reads ->
            def type

            if (meta.platform in long_platforms) {
                type = "long"
            } else if (meta.platform in short_platforms) {
                type = "short"
            } else {
                error "Unknown sequencing platform '${meta.platform}' for sample '${meta.id}'"
            }

            if (
                type == "long" &&
                reads instanceof List &&
                reads.size() > 1
            ) {
                error "Sample '${meta.id}' has platform '${meta.platform}' (long-read) but provides ${reads.size()} FASTQ files. Long-read platforms expect a single FASTQ file per sample. Check your samplesheet."
            }

            [
                meta + [
                    type: type
                ],
                reads
            ]
        }
        .set { ch_reads_mix }



    //
    // MODULE: Falco — QC on raw reads
    //
    // Only the *_fastqc_data.txt files go to MultiQC. Their names match MultiQC's
    // `fastqc/data` search pattern (`*fastqc_data.txt`), so MultiQC parses them
    // with its FastQC module and labels the section "Falco". This needs MultiQC
    // >= 1.32; 1.31 matches the filename `fastqc_data.txt` exactly and silently
    // ignores these files. The per-module pass/warn/fail statuses MultiQC shows
    // come from the data files too, so the summary files are not needed here.
    //
    FALCO (
        ch_reads_mix
    )
    ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt.collect { meta, txt -> txt })

    // MODULE: Fastp / Fastplong
    // Quality control and trimming
    //
    if (!params.skip_fastp) {

        ch_short_reads = ch_reads_mix.filter { meta, reads ->
            meta.type == 'short'
        }

        ch_long_reads = ch_reads_mix.filter { meta, reads ->
            meta.type == 'long'
        }

        FASTP (
            ch_short_reads.map { meta, reads ->
                tuple(meta, reads, [])
            },
            false,
            false,
            false
        )

        FASTPLONG (
            ch_long_reads.map { meta, reads ->
                tuple(meta, reads)
            },
            [],
            false,
            false
        )

        ch_reads_mix = FASTP.out.reads.mix(
            FASTPLONG.out.reads
        )
    }

    //
    // Group reads by sample ID + sequencing platform
    //
    ch_reads_by_sample = ch_reads_mix
        .map { meta, reads ->

            [
                [meta.id, meta.platform],
                meta,
                reads
            ]
        }
        .groupTuple(by: 0)


    //
    // Separate single-run samples from multi-run samples
    //
    ch_reads_by_sample
        .branch { key, metas, reads_list ->

            single:
                reads_list.size() == 1

                return [
                    metas[0],
                    reads_list[0]
                ]

            multiple:
                true
                // Keep each run's metadata paired with its FASTQ files,
                // then sort by the original run number before merging.
                def orderedRuns = metas.indices
                .collect { i ->
                [meta: metas[i], reads: reads_list[i]]
                }
                .sort { a, b ->
                a.meta.run <=> b.meta.run
                }

                return [
                    orderedRuns[0].meta,
                    orderedRuns.collect { it.reads }.flatten()
                    ]
                    }
                    .set { ch_reads_grouped }


    //
    // Merge samples that have multiple runs
    //
    CAT_FASTQ (
        ch_reads_grouped.multiple
    )


    //
    // Combine:
    //   - merged multi-run samples
    //   - single-run samples that do not require merging
    //
    ch_merged_reads = CAT_FASTQ.out.reads
        .mix(ch_reads_grouped.single)

    //
    // Read removal
    //
    // --dehoster selects the implementation, or turns the step off entirely.
    // With 'none' the pipeline reduces to its minimum: optional bam2fq
    // conversion, falco, and multiqc. Nothing below this point needs a host
    // index in that case, so index resolution is inside each branch — asking
    // for one when no dehosting will happen is how this used to fail.
    if (params.dehoster == 'none') {
        ch_dehosted_fastq = ch_merged_reads
        ch_dehosted_json  = channel.empty()
        ch_host_bam       = channel.empty()
        ch_host_bai       = channel.empty()
    }
    else if (params.dehoster == 'deacon') {
        // deacon is the default. It does not align: it counts how many of a
        // read's minimizers are present in a prebuilt index and drops the read
        // when enough of them are, which is where the speed over the two
        // bowtie2-based dehosters comes from.
        //
        // Two things follow from that. There is no alignment to keep, so
        // --host_bam does nothing here and hostess remains the option for
        // anyone who needs the removed reads. And the reference is a different
        // one: the index below is panhuman-1, whose masking subtracts far less
        // than the hostile index does, so the two do not remove the same reads.
        //
        // deacon's index is a single file rather than a directory of bowtie2
        // parts, so it takes a path rather than the db/name pair hostile needs.
        if (params.host_fasta) {
            def host_fasta = file(params.host_fasta, checkIfExists: true)
            DEACON_INDEX (
                channel.value([[id: host_fasta.baseName], host_fasta])
            )
            // Already a value channel, since DEACON_INDEX's input is one.
            ch_deacon_index = DEACON_INDEX.out.index.map { meta, index -> index }
        }
        else {
            // Same failure mode as hostile_db below: without this the missing
            // path surfaces as "Argument of `file()` function cannot be null",
            // which says nothing about the configs not having been loaded.
            if (!params.deacon_db) {
                error(
                    "deacon_db is not set, so there is no reference index to dehost against.\n\n" +
                    "  Pass --deacon_db, or let oamd-bio-nf-configs supply it. For the latter,\n" +
                    "  both of these must be set, and neither has a safe default:\n" +
                    "    NF_PIPELINE_REPO=oamd-bio-qqq   (otherwise it falls back to oamd-bio-ph-core,\n" +
                    "                                     loading that pipeline's configs and workDir)\n" +
                    "    NF_CONFIGS_REPO=<checkout containing conf/pipeline/oamd-bio-qqq/>\n\n" +
                    "  See run-qqq-prda.sh, which sets both.\n\n" +
                    "  --host_fasta builds an index instead, which is what the test profiles do."
                )
            }
            ch_deacon_index = channel.value(file(params.deacon_db, checkIfExists: true))
        }

        // Routed by input size into two aliases, like the aligning dehosters
        // below, but for a different reason than they are.
        //
        // Memory genuinely does not need it. deacon holds the index in memory
        // and streams reads past it, so its footprint follows the index rather
        // than the sample: on a prd-analysis run, peak RSS sat at 4.2-4.5 GB
        // for every sample from 644 KB to 3.5 GB of R1. Both tiers therefore
        // get the same memory, and only the CPU count differs.
        //
        // CPU does need it, because of how the module is built. It compresses
        // deacon's output with one bgzip per mate, and deacon and both bgzips
        // are each handed task.cpus threads, so a paired-end task asks for up
        // to three times its allocation. What it actually uses tracks the
        // sample: measured at cpus = 2, %cpu ran 107% on the 644 KB sample and
        // 451% on the 3.5 GB one, climbing monotonically in between. That is
        // the bgzip pair saturating while deacon blocks on FIFO backpressure —
        // the task is compression-bound, and compression scales with the input
        // where deacon's own filtering barely does.
        //
        // This replaces a single 2-CPU allocation whose measurement was taken
        // on 2.1M read pairs. At that size loading the index is a fixed ~6s and
        // most of the task, so extra threads multiplied a cost no thread count
        // could reduce. On the 3.5 GB sample the same load is 6s of 7m28s, and
        // the conclusion inverts with it.
        //
        // Resources are set per alias in oamd-bio-nf-configs; nothing here
        // hardcodes them, and the module itself is untouched.
        ch_deacon_sized = ch_merged_reads.branch { meta, reads ->
            large: [reads].flatten().sum { read -> read.size() } > params.hostile_large_input_bytes
            small: true
        }

        DEACON_FILTER_SMALL (
            ch_deacon_sized.small,
            ch_deacon_index
        )
        DEACON_FILTER_LARGE (
            ch_deacon_sized.large,
            ch_deacon_index
        )
        ch_dehosted_fastq = DEACON_FILTER_SMALL.out.fastq.mix(DEACON_FILTER_LARGE.out.fastq)
        ch_dehosted_json  = DEACON_FILTER_SMALL.out.json.mix(DEACON_FILTER_LARGE.out.json)
        ch_host_bam       = channel.empty()
        ch_host_bai       = channel.empty()
    }
    else {
        // A prebuilt index is the normal case. --host_fasta builds one instead,
        // which is what the test profiles use: a bowtie2 index is ~8 MB even for a
        // tiny genome, so building it is cheaper than hosting one as test data.
        if (params.host_fasta) {
            def host_fasta = file(params.host_fasta, checkIfExists: true)
            BOWTIE2_BUILD (
                channel.value([[id: host_fasta.baseName], host_fasta])
            )
            // BOWTIE2_BUILD writes <fasta.baseName>.*.bt2 into a `bowtie2` directory,
            // so the index name hostile is given must match that basename.
            // Already a value channel, since BOWTIE2_BUILD's input is one.
            ch_hostile_ref_short = BOWTIE2_BUILD.out.index.map { meta, index ->
            [meta.id, index]
            }

            // Long reads: use the original FASTA
            ch_hostile_ref_long = channel.value([
                host_fasta.name,
                host_fasta
                ])
        } else {
            // Without this, a missing hostile_db surfaces as "Argument of `file()`
            // function cannot be null", which says nothing about the actual cause:
            // the pipeline configs from oamd-bio-nf-configs were not loaded.
            if (!params.hostile_db) {
                error(
                    "hostile_db is not set, so there is no reference index to dehost against.\n\n" +
                    "  Pass --hostile_db, or let oamd-bio-nf-configs supply it. For the latter,\n" +
                    "  both of these must be set, and neither has a safe default:\n" +
                    "    NF_PIPELINE_REPO=oamd-bio-qqq   (otherwise it falls back to oamd-bio-ph-core,\n" +
                    "                                     loading that pipeline's configs and workDir)\n" +
                    "    NF_CONFIGS_REPO=<checkout containing conf/pipeline/oamd-bio-qqq/>\n\n" +
                    "  See run-qqq-prda.sh, which sets both."
                )
            }

            def hostile_db = file(params.hostile_db, checkIfExists: true)

            ch_hostile_ref_short = channel.value([
                params.hostile_index,
                hostile_db
                ])

            ch_hostile_ref_long = channel.value([
                params.hostile_index,
                hostile_db
                ])
        }

        // Hostile's peak memory plateaus rather than tracking input size: measured
        // on prd-analysis, three samples spanning 15-22 GB of R1 all peaked at
        // 43.6 GB, while samples under 100 MB peaked under 5 GB. A single
        // allocation therefore either starves the large samples or wastes tens of
        // GB on the small ones, and on a typical samplesheet most samples are
        // small.
        //
        // So route by input size to two aliases of the same module. Resources are
        // set per alias in oamd-bio-nf-configs; nothing here hardcodes them, and
        // the module itself is untouched.
        ch_sized = ch_merged_reads.branch { meta, reads ->
            large: [reads].flatten().sum { read -> read.size() } > params.hostile_large_input_bytes
            small: true
        }

        ch_small = ch_sized.small.branch { meta, reads ->
        short_reads: meta.type == 'short'
        long_reads:  meta.type == 'long'
        }

        ch_large = ch_sized.large.branch { meta, reads ->
        short_reads: meta.type == 'short'
        long_reads:  meta.type == 'long'
        }

        // hostess is the default: it removes host reads in one alignment pass and
        // can keep the removed reads as a sorted BAM, which hostile cannot do
        // without a second pass. --host_bam only has an effect with hostess.
        if (params.dehoster == 'hostess') {
            HOSTESS_SMALL_SHORT (
            ch_small.short_reads,
            ch_hostile_ref_short
            )

            HOSTESS_LARGE_SHORT (
            ch_large.short_reads,
            ch_hostile_ref_short
            )

            HOSTESS_SMALL_LONG (
            ch_small.long_reads,
            ch_hostile_ref_long
            )

            HOSTESS_LARGE_LONG (
            ch_large.long_reads,
            ch_hostile_ref_long
            )
            ch_dehosted_fastq = HOSTESS_SMALL_SHORT.out.fastq.mix(HOSTESS_LARGE_SHORT.out.fastq).mix(HOSTESS_SMALL_LONG.out.fastq).mix(HOSTESS_LARGE_LONG.out.fastq)
            ch_dehosted_json = HOSTESS_SMALL_SHORT.out.json.mix(HOSTESS_LARGE_SHORT.out.json).mix(HOSTESS_SMALL_LONG.out.json).mix(HOSTESS_LARGE_LONG.out.json)

            // Only produced when --host_bam is set; the module declares them
            // optional, so these channels are simply empty otherwise.
            ch_host_bam = HOSTESS_SMALL_SHORT.out.bam.mix(HOSTESS_LARGE_SHORT.out.bam).mix(HOSTESS_SMALL_LONG.out.bam).mix(HOSTESS_LARGE_LONG.out.bam)
            ch_host_bai = HOSTESS_SMALL_SHORT.out.bai.mix(HOSTESS_LARGE_SHORT.out.bai).mix(HOSTESS_SMALL_LONG.out.bai).mix(HOSTESS_LARGE_LONG.out.bai)

            // samtools stats is parsed by MultiQC's samtools module with no
            // configuration, so it goes straight into the report.
            //
            // flagstat is deliberately not sent. Everything MultiQC takes from it
            // duplicates stats — its "Reads mapped" was the same number as the one
            // from stats for every sample of a 49-sample run — so it bought a
            // second set of General Statistics columns and a redundant report
            // section. The file is still produced and published for anyone who
            // wants it.
            ch_multiqc_files = ch_multiqc_files.mix(
                HOSTESS_SMALL_SHORT.out.stats.mix(HOSTESS_LARGE_SHORT.out.stats).mix(HOSTESS_SMALL_LONG.out.stats).mix(HOSTESS_LARGE_LONG.out.stats).collect { meta, file -> file }.ifEmpty([])
            )
        }
        else {
            HOSTILE_CLEAN_SMALL_SHORT (
            ch_small.short_reads,
            ch_hostile_ref_short
            )

            HOSTILE_CLEAN_LARGE_SHORT (
            ch_large.short_reads,
            ch_hostile_ref_short
            )

            HOSTILE_CLEAN_SMALL_LONG (
            ch_small.long_reads,
            ch_hostile_ref_long
            )

            HOSTILE_CLEAN_LARGE_LONG (
            ch_large.long_reads,
            ch_hostile_ref_long
            )

            ch_dehosted_fastq = HOSTILE_CLEAN_SMALL_SHORT.out.fastq.mix(HOSTILE_CLEAN_LARGE_SHORT.out.fastq).mix(HOSTILE_CLEAN_SMALL_LONG.out.fastq).mix(HOSTILE_CLEAN_LARGE_LONG.out.fastq)
            ch_dehosted_json = HOSTILE_CLEAN_SMALL_SHORT.out.json.mix(HOSTILE_CLEAN_LARGE_SHORT.out.json).mix(HOSTILE_CLEAN_SMALL_LONG.out.json).mix(HOSTILE_CLEAN_LARGE_LONG.out.json)
            ch_host_bam = channel.empty()
            ch_host_bai = channel.empty()
        }
    }

    ch_multiqc_files = ch_multiqc_files.mix(ch_dehosted_json.collect { meta, json -> json }.ifEmpty([]))

    // Every dehoster writes a JSON report, and MultiQC parses it: the hostile
    // module for hostile and hostess, the deacon module for deacon. Neither
    // puts the headline number — how much of a sample was host — where a
    // reviewer looks first. The hostile module renders only a stacked bar plot
    // and contributes no General Statistics columns at all, and while the
    // deacon module does contribute some, they are its own names and scaling,
    // so the report would change shape with --dehoster.
    //
    // Reduce the reports to one MultiQC custom-content file instead, so the
    // same two columns appear on the sample's own row whichever dehoster ran.
    // The `_mqc.tsv` suffix is what makes MultiQC pick the file up, and
    // `plot_type: generalstats` in the header sends it to the top table instead
    // of a section of its own. The first column is meta.id, which is what
    // assets/multiqc_config.yml collapses the per-mate falco rows down to, so
    // these values join the existing row rather than adding another one.
    //
    // The JSON is parsed here rather than in a process because it is three fields
    // out of a few hundred bytes; a process would need a container and a task per
    // sample to do the same thing.
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_dehosted_json
            .map { meta, json -> [meta, dehostCounts(json.text.trim())] }
            // The dehoster stubs touch an empty report. A stub run has no counts
            // to put in the table, so it contributes no rows instead of rows of
            // zeroes, and with every report empty no file is written at all.
            .filter { meta, counts -> counts }
            .map { meta, counts ->
                // From the counts rather than the tool's own removed proportion,
                // which is rounded, and to keep an empty input from dividing by zero.
                def pct_host = counts.reads_in ? 100.0 * counts.reads_removed / counts.reads_in : 0.0
                [meta.id, counts.reads_out, String.format('%.4f', pct_host)].join('\t')
            }
            .collectFile(
                name: 'host_removal_mqc.tsv',
                seed: hostRemovalHeader(),
                newLine: true,
                sort: true
            )
    )

    //
    // MODULE: Fadso — cap depth of the dehosted reads (optional)
    //
    // Disabled when params.max_reads is null or 0, which is the default.
    //
    if (params.max_reads && params.dehoster == 'none') {
        // No dehoster means no read count to gate on, so every sample is routed
        // through fadso. Samples already under the cap come out unchanged, so
        // this costs a pass over the FASTQs rather than changing the result.
        FADSO (
            ch_dehosted_fastq,
            params.max_reads,
            params.downsample_seed
        )
        ch_clean_reads = FADSO.out.reads
    }
    else if (params.max_reads) {
        // Gated per sample on the dehoster's own read count, so no extra pass
        // over the FASTQs is needed. Only samples above the cap are routed to
        // FADSO; everything else bypasses the module and is merged back after.
        ch_read_counts = ch_dehosted_json
            .map { meta, json ->
                // The dehoster stubs emit an empty report. Treat that as "under
                // the cap" so -stub-run still exercises this path; a non-empty
                // but malformed report is a real error and throws.
                def counts = dehostCounts(json.text.trim())
                def reads_out = (counts ? counts.reads_out : 0L) as long
                // Every dehoster counts records rather than pairs, so for
                // paired-end input reads_out is twice the number of pairs.
                [[meta.id, meta.platform], meta.single_end ? reads_out : reads_out.intdiv(2)]
            }

        ch_gated = ch_dehosted_fastq
            .map { meta, reads -> [[meta.id, meta.platform], meta, reads] }
            .join(ch_read_counts)
            .branch { id, meta, reads, count ->
                downsample: count > params.max_reads
                    return [meta, reads]
                keep: true
                    return [meta, reads]
            }

        FADSO (
            ch_gated.downsample,
            params.max_reads,
            params.downsample_seed
        )

        ch_clean_reads = ch_gated.keep.mix(FADSO.out.reads)
    }
    else {
        ch_clean_reads = ch_dehosted_fastq
    }

    //
    // Collate and save software versions
    //
    // Every module in this pipeline emits versions to the `versions` topic
    // channel. The branch below still handles classic `versions.yml` files so a
    // module that has not been migrated yet is not silently dropped.
    //
    // MULTIQC is deliberately absent: its input depends on the versions topic
    // being resolved, so emitting to that topic would deadlock the pipeline.
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            // Strip the size-tier suffix so the two dehoster aliases report one
            // version rather than two identical ones.
            def name = process[process.lastIndexOf(':') + 1..-1].replaceAll(/_(SMALL|LARGE)_(SHORT|LONG)$/, '')
            [name, "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'qqq_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // MULTIQC takes a single tuple rather than separate channels. The pipeline's
    // own config is always passed; a user-supplied --multiqc_config is appended
    // to it rather than replacing it, so pipeline defaults still apply.
    MULTIQC (
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'multiqc'],
                files,
                params.multiqc_config
                    ? [
                        file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                        file(params.multiqc_config, checkIfExists: true),
                    ]
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    host_bam       = ch_host_bam                                      // channel: [ val(meta), path(bam) ]
    host_bai       = ch_host_bai                                      // channel: [ val(meta), path(bai) ]
    clean_reads    = ch_clean_reads                                   // channel: [ val(meta), path(fastq) ]
    multiqc_report = MULTIQC.out.report.map { meta, report -> report } // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Read counts from a dehoster's JSON report, normalised across the two
// families. hostile and hostess write a single-element array of reads_in /
// reads_out / reads_removed; deacon writes a bare object of seqs_in /
// seqs_out / seqs_removed. Both count records rather than pairs, so the
// numbers mean the same thing once the names line up.
//
// Returns null for the empty report the stubs write, so a caller can tell
// "nothing was measured" from "measured zero reads". A malformed non-empty
// report throws, which is what should happen.
//
def dehostCounts(report) {
    if (!report) {
        return null
    }
    def summary = new groovy.json.JsonSlurper().parseText(report)
    // hostile and hostess wrap theirs in an array; deacon does not.
    if (summary instanceof List) {
        summary = summary[0]
    }
    return summary.containsKey('seqs_in')
        ? [reads_in: summary.seqs_in, reads_out: summary.seqs_out, reads_removed: summary.seqs_removed]
        : [reads_in: summary.reads_in, reads_out: summary.reads_out, reads_removed: summary.reads_removed]
}

//
// Header for host_removal_mqc.tsv: the comment block MultiQC reads as custom
// content configuration, then the column names for the rows the workflow
// appends. `namespace` is what the General Statistics column configuration
// groups these two columns under, and MultiQC displays it prefixed with the
// module name, as "Custom content: Dehost".
//
def hostRemovalHeader() {
    [
        '# id: "host_removal"',
        '# plot_type: "generalstats"',
        '# namespace: "Dehost"',
        '# pconfig:',
        '#     - clean_reads:',
        '#         title: "Clean reads"',
        '#         description: "Reads remaining after host removal"',
        '#         shared_key: "read_count"',
        '#         scale: "BuGn"',
        '#     - pct_host:',
        '#         title: "% Host"',
        '#         description: "Percentage of input reads removed as host"',
        '#         min: 0',
        '#         max: 100',
        '#         suffix: "%"',
        '#         format: "{:,.2f}"',
        '#         scale: "OrRd"',
        ['Sample', 'clean_reads', 'pct_host'].join('\t'),
    ].join('\n')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

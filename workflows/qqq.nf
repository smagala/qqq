/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FALCO                  } from '../modules/amd-ph-core/falco/main'
include { FADSO                  } from '../modules/amd-ph-core/fadso/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { HOSTESS as HOSTESS_SMALL             } from '../modules/amd-ph-core/hostess/main'
include { HOSTESS as HOSTESS_LARGE             } from '../modules/amd-ph-core/hostess/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_SMALL } from '../modules/nf-core/hostile/clean/main'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_LARGE } from '../modules/nf-core/hostile/clean/main'
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
        ch_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt.collect { meta, txt -> txt })

    //
    // MODULE: Hostile clean — dehost reads
    //
    // A prebuilt index is the normal case. --hostile_fasta builds one instead,
    // which is what the test profiles use: a bowtie2 index is ~8 MB even for a
    // tiny genome, so building it is cheaper than hosting one as test data.
    if (params.hostile_fasta) {
        def hostile_fasta = file(params.hostile_fasta, checkIfExists: true)
        BOWTIE2_BUILD (
            channel.value([[id: hostile_fasta.baseName], hostile_fasta])
        )
        // BOWTIE2_BUILD writes <fasta.baseName>.*.bt2 into a `bowtie2` directory,
        // so the index name hostile is given must match that basename.
        // Already a value channel, since BOWTIE2_BUILD's input is one.
        ch_hostile_ref = BOWTIE2_BUILD.out.index.map { meta, index -> [meta.id, index] }
    }
    else {
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
        ch_hostile_ref = channel.value([params.hostile_index, file(params.hostile_db, checkIfExists: true)])
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
    ch_sized = ch_reads.branch { meta, reads ->
        large: [reads].flatten().sum { read -> read.size() } > params.hostile_large_input_bytes
        small: true
    }

    // --dehoster selects the implementation. hostess removes host reads in one
    // alignment pass and can keep the removed reads as a sorted BAM, which
    // hostile cannot do without a second pass. Both take the same inputs and
    // emit the same fastq and json channels, so everything downstream is
    // unchanged; --host_bam only has an effect with hostess.
    if (params.dehoster == 'hostess') {
        HOSTESS_SMALL (
            ch_sized.small,
            ch_hostile_ref
        )
        HOSTESS_LARGE (
            ch_sized.large,
            ch_hostile_ref
        )
        ch_hostile_fastq = HOSTESS_SMALL.out.fastq.mix(HOSTESS_LARGE.out.fastq)
        ch_hostile_json = HOSTESS_SMALL.out.json.mix(HOSTESS_LARGE.out.json)

        // Only produced when --host_bam is set; the module declares them
        // optional, so these channels are simply empty otherwise.
        ch_host_bam = HOSTESS_SMALL.out.bam.mix(HOSTESS_LARGE.out.bam)
        ch_host_bai = HOSTESS_SMALL.out.bai.mix(HOSTESS_LARGE.out.bai)

        // samtools stats and flagstat are parsed by MultiQC's samtools module
        // with no configuration, so they go straight into the report.
        ch_multiqc_files = ch_multiqc_files.mix(
            HOSTESS_SMALL.out.flagstat.mix(HOSTESS_LARGE.out.flagstat).collect { meta, file -> file }.ifEmpty([])
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            HOSTESS_SMALL.out.stats.mix(HOSTESS_LARGE.out.stats).collect { meta, file -> file }.ifEmpty([])
        )
    }
    else {
        HOSTILE_CLEAN_SMALL (
            ch_sized.small,
            ch_hostile_ref
        )
        HOSTILE_CLEAN_LARGE (
            ch_sized.large,
            ch_hostile_ref
        )
        ch_hostile_fastq = HOSTILE_CLEAN_SMALL.out.fastq.mix(HOSTILE_CLEAN_LARGE.out.fastq)
        ch_hostile_json = HOSTILE_CLEAN_SMALL.out.json.mix(HOSTILE_CLEAN_LARGE.out.json)
        ch_host_bam = channel.empty()
        ch_host_bai = channel.empty()
    }

    ch_multiqc_files = ch_multiqc_files.mix(ch_hostile_json.collect { meta, json -> json })

    //
    // MODULE: Fadso — cap depth of the dehosted reads (optional)
    //
    // Gated per sample on HOSTILE_CLEAN's own read count, so no extra pass over
    // the FASTQs is needed. Only samples above the cap are routed to FADSO;
    // everything else bypasses the module entirely and is merged back after.
    // Disabled when params.max_reads is null or 0.
    //
    if (params.max_reads) {
        ch_read_counts = ch_hostile_json
            .map { meta, json ->
                // HOSTILE_CLEAN's stub emits an empty report. Treat that as
                // "under the cap" so -stub-run still exercises this path; a
                // non-empty but malformed report is a real error and throws.
                def report = json.text.trim()
                def reads_out = report
                    ? new groovy.json.JsonSlurper().parseText(report)[0].reads_out as long
                    : 0L
                // Hostile counts alignment records (`samtools view -F 2304 -c`),
                // so for paired-end input reads_out is twice the number of pairs.
                [meta.id, meta.single_end ? reads_out : reads_out.intdiv(2)]
            }

        ch_gated = ch_hostile_fastq
            .map { meta, reads -> [meta.id, meta, reads] }
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
        ch_clean_reads = ch_hostile_fastq
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
            // Strip the size-tier suffix so the two HOSTILE_CLEAN aliases report
            // one hostile version rather than two identical ones.
            def name = process[process.lastIndexOf(':') + 1..-1].replaceAll(/_(SMALL|LARGE)$/, '')
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
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

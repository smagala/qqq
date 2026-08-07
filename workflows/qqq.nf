/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FALCO                  } from '../modules/amd-ph-core/falco/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { HOSTILE_CLEAN          } from '../modules/nf-core/hostile/clean/main'
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
        ch_hostile_ref = channel.value([params.hostile_index, file(params.hostile_db, checkIfExists: true)])
    }

    HOSTILE_CLEAN (
        ch_reads,
        ch_hostile_ref
    )
    ch_multiqc_files = ch_multiqc_files.mix(HOSTILE_CLEAN.out.json.collect { meta, json -> json })

    ch_clean_reads = HOSTILE_CLEAN.out.fastq

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
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
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
    clean_reads    = ch_clean_reads                                   // channel: [ val(meta), path(fastq) ]
    multiqc_report = MULTIQC.out.report.map { meta, report -> report } // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

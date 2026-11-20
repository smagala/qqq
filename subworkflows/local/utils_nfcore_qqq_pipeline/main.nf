//
// Subworkflow with functionality specific to the cdcent/oamd-bio-qqq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null    // cli_typecast: nf-schema's default, which casts command-line values using the schema. On
                // Nextflow 26.04 they arrive as strings; `false` made every boolean or integer flag fail validation.
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    // ph-core allows only `id` and `single_end` in meta, so the sequencing
    // platform and the run number are not meta fields. Both are still needed:
    // the platform routes short- against long-read tooling, and the run number
    // orders a sample's runs before they are merged.
    //
    // They travel in a `route` map that rides beside the meta in the channel and
    // is dropped before any process call -- see workflows/qqq.nf. A value that
    // never enters a module cannot leak into one, which is a stronger guarantee
    // than a meta key that merely happens to be unused downstream.
    //
    // `route.merge_id` is the identity a sample has once its runs are merged.
    // It is carried rather than recovered by stripping `_run<n>` off `meta.id`:
    // parsing identity back out of an id is how ids stop being opaque.
    def sheet = samplesheetToList(input, "${projectDir}/assets/schema_input.json")

    // merge_id is DERIVED from two fields, so two rows can derive the same
    // string without being the same sample: the platform enum holds both
    // `nanopore` and `oxford_nanopore`, which makes sample `x_oxford` on
    // nanopore and sample `x` on oxford_nanopore both "x_oxford_nanopore".
    // Grouping on the derived string would merge two different samples' reads
    // into one FASTQ, and publishing on it would have them overwrite each
    // other's outputs -- silently, in both cases.
    //
    // So nothing groups on the derived string. The grouping key below is the
    // structured pair, and the string is checked for uniqueness here, over the
    // whole sheet, before any of it reaches a channel: a collision is a
    // property of the samplesheet, knowable before the first task is submitted,
    // and better failed on than carried into the DAG.
    def merge_id_source = [:]
    sheet.each { entry ->
        def row      = entry[0]
        def platform = normalisePlatform(row)
        def merge_id = "${row.id}_${platform}"
        def source   = [ row.id, platform ]
        if (merge_id_source.containsKey(merge_id) && merge_id_source[merge_id] != source) {
            def other = merge_id_source[merge_id]
            error(
                "Samplesheet rows (id '${other[0]}', platform '${other[1]}') and " +
                "(id '${source[0]}', platform '${source[1]}') are different samples but both " +
                "merge to the identity '${merge_id}', so their reads would be merged together " +
                "and their outputs would collide. Rename one of the samples."
            )
        }
        merge_id_source[merge_id] = source
    }

    channel
        .fromList(sheet)
        .map {
            // `row` is the samplesheet row as the schema defines it, not pipeline
            // meta. Named apart from `meta` so that reading row.platform here is
            // not mistaken for a meta key, which is what this change removes.
            row, fastq_1, fastq_2 ->
            def platform = normalisePlatform(row)
            def merge_id = "${row.id}_${platform}"
            def fastqs   = fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ]

            // meta is built explicitly rather than spread from the samplesheet
            // row. Spreading is how `platform` became a meta key in the first
            // place: every column the schema gains arrives in meta for free.
            return [
                [ row.id, platform ],
                merge_id,
                [ id: merge_id, single_end: !fastq_2 ],
                fastqs
            ]
        }
        .groupTuple(by: 0)
        .map { source, merge_ids, metas, fastqs ->
            // Both are constant within a group: the group key is exactly what
            // they are derived from, so any element identifies the group.
            def platform = source[1]
            def merge_id = merge_ids[0]
            validateInputSamplesheet([merge_id, metas, fastqs])
            [merge_id, platform, metas, fastqs]
        }
        .flatMap { merge_id, platform, metas, fastqs ->
            // `runs` is the number of runs this sample will be merged from.
            // It is emitted here because here is the only place that knows it
            // without waiting: the sheet is a finite list, so the group is
            // already complete. Downstream, qqq regroups runs into samples and
            // needs the group size to do it incrementally -- see the groupKey
            // in workflows/qqq.nf. Recovering it later would mean counting the
            // runs as they arrive, which is the thing being avoided.
            def runs = metas.size()
            metas.indices.collect { i ->
                def run = i + 1
                [
                    [ id: "${merge_id}_run${run}", single_end: metas[i].single_end ],
                    [ platform: platform, run: run, runs: runs, merge_id: merge_id ],
                    fastqs[i]
                ]
            }
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
}

//
// Normalise the samplesheet's platform column.
//
// A function because the value is read twice -- once to check that the merged
// identities it produces are unique, once to build them -- and two copies of
// the normalisation would be two chances for the check to normalise differently
// from the thing it is checking.
//
def normalisePlatform(row) {
    return row.platform?.toLowerCase()?.trim() ?: "illumina"
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (_merge_id, metas, fastqs) = input

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    def all_fastqs = fastqs.flatten()

    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
        }

    if (all_fastqs.size() != all_fastqs.unique().size()) {
        error( "Please check input samplesheet -> Duplicate FASTQ files detected for sample: ${metas[0].id}")
        }

    // A sample may be split over several files (lanes, re-sequencing runs) and they are merged.
    // But the same file stored twice -- e.g. a BaseSpace dataset copy (`.../ds.<id>/<lane file>`)
    // listed alongside the original -- would merge every read in it twice. Two files with the same
    // name and the same size are taken to be one file stored twice. The same name with a different
    // size is allowed: a re-sequencing run of one library reuses the instrument's file names.
    def copies = all_fastqs
        .groupBy { fq -> [ fq.name, fq.size() ] }
        .findAll { _key, files -> files.size() > 1 }
    if (copies) {
        def listed = copies.collect { key, files -> "  ${key[0]} (${key[1]} bytes):\n" + files.collect { f -> "    ${f.toUriString()}" }.join("\n") }.join("\n")
        error( "Please check input samplesheet -> sample ${metas[0].id} lists the same file more than once (same name and size, different paths); its reads would be merged twice. Keep one copy of each:\n${listed}")
        }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
//
// Only the tools a given run actually used are named. The dehoster is chosen
// at runtime, so citing all three would credit two tools that never ran, and
// Bowtie 2 only appears when an aligning dehoster built its index from a FASTA.
//
def toolCitationText() {
    def aligning_dehoster = params.dehoster in ['hostile', 'hostess']
    def citation_text = [
            "Tools used in the workflow included:",
            "Falco (de Sena Brandine and Smith 2019),",
            !params.skip_fastp ? "fastp (Chen et al. 2018) and fastplong (Chen 2023), as applicable," : "",
            aligning_dehoster && params.host_fasta ? "Bowtie 2 (Langmead and Salzberg 2012)," : "",
            params.dehoster == 'deacon' ? "Deacon (Constantinides et al. 2025)," : "",
            aligning_dehoster ? "Hostile (Constantinides et al. 2023)," : "",
            params.max_reads ? "Rasusa (Hall 2022)," : "",
            "and MultiQC (Ewels et al. 2016)."
        ].findAll { entry -> entry }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def aligning_dehoster = params.dehoster in ['hostile', 'hostess']
    def reference_text = [
            "<li>de Sena Brandine, G., & Smith, A. D. (2019). Falco: high-speed FastQC emulation for quality control of sequencing data. F1000Research, 8, 1874. doi: 10.12688/f1000research.21142.2</li>",
            !params.skip_fastp ? "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. doi: 10.1093/bioinformatics/bty560</li>" : "",
            !params.skip_fastp ? "<li>Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta, 2, e107. doi: 10.1002/imt2.107</li>" : "",
            aligning_dehoster && params.host_fasta ? "<li>Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357-359. doi: 10.1038/nmeth.1923</li>" : "",
            params.dehoster == 'deacon' ? "<li>Constantinides, B., Lees, J., & Crook, D. W. (2025). Deacon: fast sequence filtering and contaminant depletion. bioRxiv. doi: 10.1101/2025.06.09.658732</li>" : "",
            aligning_dehoster ? "<li>Constantinides, B., Hunt, M., & Crook, D. W. (2023). Hostile: accurate decontamination of microbial host sequences. Bioinformatics, 39(12), btad728. doi: 10.1093/bioinformatics/btad728</li>" : "",
            params.max_reads ? "<li>Hall, M. B. (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941. doi: 10.21105/joss.03941</li>" : "",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
        ].findAll { entry -> entry }.join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

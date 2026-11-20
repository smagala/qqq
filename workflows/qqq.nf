/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FALCO                  } from '../modules/amd-ph-core/falco/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTPLONG              } from '../modules/nf-core/fastplong/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { RASUSA_READS           } from '../modules/amd-ph-core/rasusa/reads/main'
include { BBMAP_BBDUK            } from '../modules/amd-ph-core/bbmap/bbduk/main'
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
include { QQQ_RUN_REPORT         } from '../modules/local/qqq_run_report/main'
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

    // `route` carries what ph-core does not allow in meta: the platform, the run
    // number and the merged-sample identity. It is an ordinary channel element
    // that is dropped at every process call, so it cannot reach a module -- a
    // stronger guarantee than a meta key that merely happens to go unread.
    // `type` is derived here rather than in the parser because the platform
    // lists it is derived from live in this file.
    ch_reads
        .map { meta, route, reads ->
            def type

            if (route.platform in long_platforms) {
                type = "long"
            } else if (route.platform in short_platforms) {
                type = "short"
            } else {
                error "Unknown sequencing platform '${route.platform}' for sample '${meta.id}'"
            }

            if (
                type == "long" &&
                reads instanceof List &&
                reads.size() > 1
            ) {
                error "Sample '${meta.id}' has platform '${route.platform}' (long-read) but provides ${reads.size()} FASTQ files. Long-read platforms expect a single FASTQ file per sample. Check your samplesheet."
            }

            [
                meta,
                route + [ type: type ],
                reads
            ]
        }
        .set { ch_reads_mix }

    // Kept as it is here, before fastp replaces the reads and before runs are
    // merged, for the run report: it is the only point where each run's
    // samplesheet inputs and the sample it will be merged into sit together.
    ch_input_runs = ch_reads_mix

    // Records for QQQ_RUN_REPORT, filled in as each step is wired below. Each is
    // [ record map, [ files the report task must stage ] ]; see the report
    // section near the end of this workflow for what the records mean.
    ch_report_records = ch_input_runs.map { meta, route, reads ->
        [
            [
                kind: 'input_run', sample: route.merge_id, run: meta.id, run_number: route.run,
                single_end: meta.single_end, platform: route.platform, read_type: route.type,
                files: [reads].flatten().collect { f -> f.toUriString() },
            ],
            [],
        ]
    }
    // The lowest-priority answer to "which step last wrote this sample's reads":
    // none, so they are the samplesheet inputs. Only a single-run sample with
    // every step turned off ends up here. Nextflow publishes only files a task
    // wrote, never pipeline inputs, so reads/ has no copy of these; the report
    // gives the samplesheet paths instead.
    ch_report_records = ch_report_records.mix(
        ch_input_runs.map { _meta, route, reads ->
            [
                [
                    kind: 'reads', sample: route.merge_id, step: 'input', priority: 0, published: false,
                    source_files: [reads].flatten().collect { f -> f.toUriString() },
                ],
                [],
            ]
        }
    )
    // run id -> the sample it is merged into, for re-keying per-run fastp output.
    ch_run_to_sample = ch_input_runs.map { meta, route, _reads -> [ meta.id, route.merge_id ] }


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
        ch_reads_mix.map { meta, _route, reads -> [ meta, reads ] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt.collect { _meta, txt -> txt })

    // MODULE: Fastp / Fastplong
    // Quality control and trimming
    //
    if (!params.skip_fastp) {

        ch_short_reads = ch_reads_mix.filter { _meta, route, _reads ->
            route.type == 'short'
        }

        ch_long_reads = ch_reads_mix.filter { _meta, route, _reads ->
            route.type == 'long'
        }

        // --save_failed_reads hands fastp `--failed_out` (plus `--unpaired1/2` for
        // paired-end, which catch a mate that passed when its partner did not --
        // the pair is dropped together, so those reads are removed too). They are
        // published beside the trimmed reads. On by default because the filter
        // drops a few percent at most, so keeping them is cheap, and without them
        // an over-strict filter is visible only as a count.
        //
        // They are not the whole of what paired-end fastp removes. A pair in which
        // BOTH mates fail is counted as removed and written nowhere -- measured on
        // fastp 1.3.6, 32 of 373 removed pairs on the test sample, with or without
        // --unpaired1/2. Single-end fastp and fastplong write every failed read.
        // The run report marks the paired-end entries incomplete for this reason.
        FASTP (
            ch_short_reads.map { meta, _route, reads ->
                tuple(meta, reads, [])
            },
            false,
            params.save_failed_reads,
            false
        )

        FASTPLONG (
            ch_long_reads.map { meta, _route, reads ->
                tuple(meta, reads)
            },
            [],
            false,
            params.save_failed_reads
        )

        // For the run report. Per run, re-keyed to the sample the run merges
        // into. `fastp/` and `fastplong/` are where the default publishDir rule
        // in conf/modules.config puts these processes' output.
        ch_report_records = ch_report_records
            .mix(
                FASTP.out.json.map { meta, json -> [ meta.id, 'fastp', json ] }
                    .mix(FASTPLONG.out.json.map { meta, json -> [ meta.id, 'fastplong', json ] })
                    .combine(ch_run_to_sample, by: 0)
                    .map { run, tool, json, sample ->
                        [ [ kind: 'fastp_report', sample: sample, run: run, tool: tool, file: json.name ], [ json ] ]
                    }
            )
            .mix(
                FASTP.out.reads_fail.map { meta, fail -> [ meta.id, 'fastp', fail ] }
                    .mix(FASTPLONG.out.reads_fail.map { meta, fail -> [ meta.id, 'fastplong', fail ] })
                    .combine(ch_run_to_sample, by: 0)
                    .map { run, dir, fail, sample ->
                        [ [ kind: 'fastp_failed', sample: sample, run: run, files: reportPublished(dir, fail) ], [] ]
                    }
            )
            .mix(
                FASTP.out.reads.map { meta, reads -> [ meta.id, 'fastp', reads ] }
                    .mix(FASTPLONG.out.reads.map { meta, reads -> [ meta.id, 'fastplong', reads ] })
                    .combine(ch_run_to_sample, by: 0)
                    .map { _run, _dir, reads, sample ->
                        [ reportReads(sample, 'read_filter', 1, reads), [] ]
                    }
            )

        // fastp drops `route`, so it is rejoined on meta.id. The join is keyed
        // rather than re-derived: at this point each id is still per-run and
        // unique, and fastp emits one record per input.
        ch_route_by_id = ch_reads_mix.map { meta, route, _reads -> [ meta.id, route ] }

        ch_reads_mix = FASTP.out.reads
            .mix(FASTPLONG.out.reads)
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_route_by_id)
            .map { _id, meta, reads, route -> [ meta, route, reads ] }
    }

    //
    // Group reads by sample ID + sequencing platform
    //
    // route.merge_id is already "<sample>_<platform>", so it is the grouping key
    // on its own -- no composite key, and nothing read out of meta.
    //
    // The key is wrapped in groupKey with route.runs, the number of runs this
    // sample is merged from, because an unsized groupTuple cannot emit any
    // group until its upstream channel is complete. Unsized, that makes the
    // slowest FASTP in the run a barrier in front of every downstream task for
    // every other sample: on a 303-sample run, 302 FASTP tasks finished and
    // zero BBMAP_BBDUK tasks had been submitted, because one 62 GB library
    // among ~14 GB peers was still trimming. The run executes as two serial
    // phases rather than as a pipeline, and Batch usage is correspondingly
    // spiky. With the size known per group, each sample flows through as soon
    // as its own runs are done.
    //
    // groupKey rather than `size:`, which takes one constant for every group
    // and would be wrong for any sheet mixing single-run and multi-run samples.
    ch_reads_by_sample = ch_reads_mix
        .map { meta, route, reads ->

            [
                groupKey(route.merge_id, route.runs),
                meta,
                route,
                reads
            ]
        }
        .groupTuple(by: 0)
        .map { _key, metas, routes, reads_list ->
            // Unwrap the GroupKey back to the plain merge_id string. It is read
            // off route rather than from key.toString() because meta.id is
            // built from it below and then joined on downstream: an id has to
            // be a String, not a GroupKey that prints like one.
            [ routes[0].merge_id, metas, routes, reads_list ]
        }


    //
    // Separate single-run samples from multi-run samples
    //
    // Both branches re-key meta to the merged identity, dropping the `_run<n>`
    // suffix that falco and fastp were named by. A single-run sample needs this
    // as much as a merged one: everything downstream of here is named for the
    // sample, not for the run it happened to arrive in.
    //
    // `run` is deliberately not carried past this point. It describes an input
    // file, and once the runs are merged there is no single run to describe.
    ch_reads_by_sample
        .branch { merge_id, metas, routes, reads_list ->

            single:
                reads_list.size() == 1

                return [
                    [ id: merge_id, single_end: metas[0].single_end ],
                    [ platform: routes[0].platform, type: routes[0].type, merge_id: merge_id ],
                    reads_list[0]
                ]

            multiple:
                true
                // Keep each run's route paired with its FASTQ files, then sort
                // by the original run number before merging, so the merged file
                // is in samplesheet order rather than in channel-arrival order.
                def orderedRuns = routes.indices
                    .collect { i ->
                        [ route: routes[i], reads: reads_list[i] ]
                    }
                    .sort { a, b ->
                        a.route.run <=> b.route.run
                    }

                return [
                    [ id: merge_id, single_end: metas[0].single_end ],
                    [ platform: routes[0].platform, type: routes[0].type, merge_id: merge_id ],
                    orderedRuns.collect { run -> run.reads }.flatten()
                ]
        }
        .set { ch_reads_grouped }


    //
    // Merge samples that have multiple runs
    //
    CAT_FASTQ (
        ch_reads_grouped.multiple.map { meta, _route, reads -> [ meta, reads ] }
    )
    ch_report_records = ch_report_records.mix(
        CAT_FASTQ.out.reads.map { meta, reads -> [ reportReads(meta.id, 'merge', 2, reads), [] ] }
    )


    //
    // Combine:
    //   - merged multi-run samples
    //   - single-run samples that do not require merging
    //
    // CAT_FASTQ drops `route`, so it is rejoined on meta.id, which is now the
    // merged identity and unique per sample.
    ch_merged_reads = CAT_FASTQ.out.reads
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .join(ch_reads_grouped.multiple.map { meta, route, _reads -> [ meta.id, route ] })
        .map { _id, meta, reads, route -> [ meta, route, reads ] }
        .mix(ch_reads_grouped.single)

    //
    // MODULE: BBDuk — ribosomal RNA depletion (optional)
    //
    // Off unless --rrna_ref names a panel: there is no sensible built-in
    // reference, and silently depleting nothing would be worse than not running.
    //
    // It runs BEFORE dehosting, and before downsampling.
    //
    // Before dehosting, because rRNA is the most conserved sequence there is, so
    // a host index takes more than host rRNA. deacon removes a read on shared
    // minimizers, and conserved 18S minimizers are human ones: on simulated
    // reads it removed 51-53% of fungal 18S (Candida albicans, Aspergillus
    // fumigatus), 10% of Plasmodium falciparum 18S and 9% of nematode 18S, and
    // keeps none of what it removes. Run first, bbduk takes all of it, from
    // every domain, and keeps it as the matched stream, so non-human eukaryotic
    // rRNA is characterisable instead of silently gone. It also makes the two
    // report columns mean what they say: "rRNA" is rRNA of any origin, human
    // included, and "host" is host reads that are not rRNA.
    //
    // The price is that bbduk now reads the host reads too. On the 198-sample
    // arbovirus mNGS run (2026-11-08) bbduk's realtime fit 47 s + 10 s per Gbase
    // of input (r = 0.89), the 47 s being the reference load, which does not
    // move. Moving it ahead of deacon adds an estimated 0.38 CPU-hours to the
    // run's 2.98, a median of 0.2 min for samples over half host, and 2.0 min
    // for the worst, a 13.6 Gbase library that was 79% host. Peak RSS follows
    // the reference, not the reads.
    //
    // A consequence for the matched stream: it now carries host rRNA, including
    // mitochondrial 12S/16S, and the mate of any read that hit -- bbduk moves
    // the pair. Treat bbduk/*.matched* as containing human sequence.
    //
    // Before downsampling, because after it the cap would already have been
    // spent on ribosomes, which on an RNA library is most of the sample -- 94%
    // of the wastewater library this was built against -- so a 1M-read cap would
    // leave ~60k reads of anything else.
    //
    // k-mer matching rather than alignment, because on a panel of tens of
    // thousands of near-identical rRNA records a minimizer aligner discards
    // exactly the seeds that identify the right record: measured at 24.5% of a
    // library against blastn's 68.4% at the identical identity threshold, where
    // bbduk k=25 takes 94.1%. The module's meta.yml has the numbers.
    //
    // The reference to point this at is a versioned, checksummed build-time
    // asset, not a panel assembled per run; see the rrna_ref help text for the
    // one AMDP publishes and what it removes. Its README is the authority on
    // what the per-record outputs may and may not be read as.
    //
    // The matched reads are kept rather than dropped. They are the input to rRNA
    // characterisation, and they are also the only way an over-aggressive panel
    // shows up in the output instead of as a quietly smaller library.
    //
    // bbduk takes plain [meta, reads], and the hostile/hostess routing below
    // needs `route` (short vs long), so it is rejoined on meta.id afterwards.
    //
    if (params.rrna_ref) {
        BBMAP_BBDUK (
            ch_merged_reads.map { meta, _route, reads -> [ meta, reads ] },
            file(params.rrna_ref, checkIfExists: true)
        )
        ch_predehost_reads = BBMAP_BBDUK.out.reads
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_merged_reads.map { meta, route, _reads -> [ meta.id, route ] })
            .map { _id, meta, reads, route -> [ meta, route, reads ] }
        ch_rrna_reads     = BBMAP_BBDUK.out.matched
        // MultiQC reads this with its `bbmap/stats` pattern, which content-matches
        // the file's own header rather than a filename or a tool banner. bbduk's
        // stderr log is deliberately not wired in: MultiQC's `bbduk` pattern looks
        // for "Executing jgi.BBDuk", which BBTools 39.76 renamed to
        // `bbduk.BBDukS`, and which the JVM puts on line three where the pattern
        // allows two. It matches neither version.
        // The TRUNCATED stats copy, not the full one. MultiQC writes every row
        // of what it parses into its own outputs, and the full file has one row
        // per matched reference record -- 93,601 of them against SILVA SSU+LSU
        // NR99, which costs 23 MB of multiqc_data.json for a single sample. The
        // truncated copy parses identically (the `bbmap/stats` pattern matches
        // on the header) and gives the same percent-filtered figure from 38 KB.
        // The full file and the per-record rpkm stay in the results.
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.stats_mqc.collect { _meta, txt -> txt })

        // For the run report: the stderr log, whose `Total Removed:`/`Result:` lines
        // are the exact record counts, and the complete stats file for #Total and
        // #Matched. #Matched alone undercounts paired-end removal (one matching mate
        // removes the pair). `bbduk/` is the explicit publishDir for BBMAP_BBDUK in
        // conf/modules.config.
        ch_report_records = ch_report_records
            .mix(BBMAP_BBDUK.out.stats.join(BBMAP_BBDUK.out.log).map { meta, txt, log ->
                [ [ kind: 'rrna_report', sample: meta.id, file: txt.name, log: log.name ], [ txt, log ] ]
            })
            .mix(BBMAP_BBDUK.out.matched.map { meta, reads ->
                [ [ kind: 'rrna_matched', sample: meta.id, files: reportPublished('bbduk', reads) ], [] ]
            })
            .mix(BBMAP_BBDUK.out.reads.map { meta, reads ->
                [ reportReads(meta.id, 'rrna_depletion', 3, reads), [] ]
            })
    }
    else {
        ch_predehost_reads = ch_merged_reads
        ch_rrna_reads      = channel.empty()
    }


    //
    // Read removal
    //
    // --dehoster selects the implementation, or turns the step off entirely.
    // With 'none' the pipeline reduces to its minimum: optional bam2fq
    // conversion, falco, and multiqc. Nothing below this point needs a host
    // index in that case, so index resolution is inside each branch — asking
    // for one when no dehosting will happen is how this used to fail.
    if (params.dehoster == 'none') {
        ch_dehosted_fastq = ch_predehost_reads.map { meta, _route, reads -> [ meta, reads ] }
        ch_dehosted_json  = channel.empty()
        ch_host_bam       = channel.empty()
        ch_host_bai       = channel.empty()
        ch_host_clean     = channel.empty()
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
            ch_deacon_index = DEACON_INDEX.out.index.map { _meta, index -> index }
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
        // deacon does not care short from long, so `route` is dropped on the way
        // in and the branch is on input size alone.
        ch_deacon_sized = ch_predehost_reads
            .map { meta, _route, reads -> [ meta, reads ] }
            .branch { _meta, reads ->
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
        ch_host_clean     = ch_dehosted_fastq
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
        // Sized first, then split short from long. `route` is carried through
        // both branches and dropped on the way into each alias, so the four
        // aliases below receive the plain [meta, reads] every module expects.
        ch_sized = ch_predehost_reads.branch { _meta, _route, reads ->
            large: [reads].flatten().sum { read -> read.size() } > params.hostile_large_input_bytes
            small: true
        }

        ch_small = ch_sized.small
            .branch { _meta, route, _reads ->
                short_reads: route.type == 'short'
                long_reads:  route.type == 'long'
            }

        ch_large = ch_sized.large
            .branch { _meta, route, _reads ->
                short_reads: route.type == 'short'
                long_reads:  route.type == 'long'
            }

        // Drop `route` once the routing it exists for has happened.
        ch_small_short = ch_small.short_reads.map { meta, _route, reads -> [ meta, reads ] }
        ch_small_long  = ch_small.long_reads.map  { meta, _route, reads -> [ meta, reads ] }
        ch_large_short = ch_large.short_reads.map { meta, _route, reads -> [ meta, reads ] }
        ch_large_long  = ch_large.long_reads.map  { meta, _route, reads -> [ meta, reads ] }

        // hostess is the default: it removes host reads in one alignment pass and
        // can keep the removed reads as a sorted BAM, which hostile cannot do
        // without a second pass. --host_bam only has an effect with hostess.
        if (params.dehoster == 'hostess') {
            // The aligner is passed per alias rather than derived inside the
            // module: bowtie2 takes an index prefix, minimap2 a single FASTA or
            // .mmi, and these four aliases already separate long from short with
            // their own reference channels. The module used to read meta.type for
            // this, which ph-core forbids -- only meta.id and meta.single_end are
            // allowed meta fields.
            HOSTESS_SMALL_SHORT (
            ch_small_short,
            ch_hostile_ref_short,
            'bowtie2'
            )

            HOSTESS_LARGE_SHORT (
            ch_large_short,
            ch_hostile_ref_short,
            'bowtie2'
            )

            HOSTESS_SMALL_LONG (
            ch_small_long,
            ch_hostile_ref_long,
            'minimap2'
            )

            HOSTESS_LARGE_LONG (
            ch_large_long,
            ch_hostile_ref_long,
            'minimap2'
            )
            ch_dehosted_fastq = HOSTESS_SMALL_SHORT.out.fastq.mix(HOSTESS_LARGE_SHORT.out.fastq).mix(HOSTESS_SMALL_LONG.out.fastq).mix(HOSTESS_LARGE_LONG.out.fastq)
            ch_dehosted_json = HOSTESS_SMALL_SHORT.out.json.mix(HOSTESS_LARGE_SHORT.out.json).mix(HOSTESS_SMALL_LONG.out.json).mix(HOSTESS_LARGE_LONG.out.json)

            // Only produced when --host_bam is set; the module declares them
            // optional, so these channels are simply empty otherwise.
            ch_host_bam = HOSTESS_SMALL_SHORT.out.bam.mix(HOSTESS_LARGE_SHORT.out.bam).mix(HOSTESS_SMALL_LONG.out.bam).mix(HOSTESS_LARGE_LONG.out.bam)
            ch_host_bai = HOSTESS_SMALL_SHORT.out.bai.mix(HOSTESS_LARGE_SHORT.out.bai).mix(HOSTESS_SMALL_LONG.out.bai).mix(HOSTESS_LARGE_LONG.out.bai)
            ch_host_clean = ch_dehosted_fastq

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
                HOSTESS_SMALL_SHORT.out.stats.mix(HOSTESS_LARGE_SHORT.out.stats).mix(HOSTESS_SMALL_LONG.out.stats).mix(HOSTESS_LARGE_LONG.out.stats).collect { _meta, file -> file }.ifEmpty([])
            )
        }
        else {
            // As with hostess, the aligner is an input rather than a meta field.
            // The patched module used to read meta.type to pick between a
            // bowtie2 index prefix and a minimap2 FASTA/.mmi, and to pick the
            // hostile cache directory; it now takes the aligner and emits
            // --aligner itself, so it is named in one place.
            HOSTILE_CLEAN_SMALL_SHORT (
            ch_small_short,
            ch_hostile_ref_short,
            'bowtie2'
            )

            HOSTILE_CLEAN_LARGE_SHORT (
            ch_large_short,
            ch_hostile_ref_short,
            'bowtie2'
            )

            HOSTILE_CLEAN_SMALL_LONG (
            ch_small_long,
            ch_hostile_ref_long,
            'minimap2'
            )

            HOSTILE_CLEAN_LARGE_LONG (
            ch_large_long,
            ch_hostile_ref_long,
            'minimap2'
            )

            ch_dehosted_fastq = HOSTILE_CLEAN_SMALL_SHORT.out.fastq.mix(HOSTILE_CLEAN_LARGE_SHORT.out.fastq).mix(HOSTILE_CLEAN_SMALL_LONG.out.fastq).mix(HOSTILE_CLEAN_LARGE_LONG.out.fastq)
            ch_dehosted_json = HOSTILE_CLEAN_SMALL_SHORT.out.json.mix(HOSTILE_CLEAN_LARGE_SHORT.out.json).mix(HOSTILE_CLEAN_SMALL_LONG.out.json).mix(HOSTILE_CLEAN_LARGE_LONG.out.json)
            ch_host_bam = channel.empty()
            ch_host_bai = channel.empty()
            ch_host_clean = ch_dehosted_fastq
        }
    }

    ch_multiqc_files = ch_multiqc_files.mix(ch_dehosted_json.collect { _meta, json -> json }.ifEmpty([]))

    // For the run report. The dehoster's directory is its own name: that is what
    // the default publishDir rule in conf/modules.config derives from each
    // alias (DEACON_FILTER_SMALL -> deacon, HOSTILE_CLEAN_LARGE_LONG -> hostile).
    // The dehosted reads themselves go to reads/, not here. Host reads are published
    // only as hostess's BAM; deacon and hostile keep none, deliberately.
    ch_report_records = ch_report_records
        .mix(ch_dehosted_json.map { meta, json ->
            [ [ kind: 'host_report', sample: meta.id, tool: params.dehoster, file: json.name ], [ json ] ]
        })
        .mix(ch_host_clean.map { meta, reads ->
            [ reportReads(meta.id, 'host_removal', 4, reads), [] ]
        })
        .mix(ch_host_bam.mix(ch_host_bai).map { meta, bam ->
            [ [ kind: 'host_bam', sample: meta.id, files: reportPublished('hostess', bam) ], [] ]
        })

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
            .filter { _meta, counts -> counts }
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
    // MODULE: rasusa — cap depth of the dehosted reads (optional)
    //
    // Disabled when params.max_reads is null or 0, which is the default.
    //
    if (params.max_reads && params.dehoster == 'none') {
        // No dehoster means no read count to gate on, so every sample is routed
        // through rasusa. Samples already under the cap come out whole, so this
        // costs a pass over the FASTQs rather than changing the result.
        RASUSA_READS (
            ch_dehosted_fastq,
            params.max_reads,
            params.downsample_seed
        )
        ch_clean_reads = RASUSA_READS.out.reads
    }
    else if (params.max_reads) {
        // Gated per sample on the dehoster's own read count, so no extra pass
        // over the FASTQs is needed. Only samples above the cap are routed to
        // rasusa; everything else bypasses the module and is merged back after.
        // The gate is an optimisation rather than a correctness requirement --
        // an input under the cap comes out of rasusa whole -- so a count that
        // has gone stale costs a pass over the FASTQs and nothing else. With
        // rRNA depletion ahead of dehosting, the dehoster is the last step to
        // touch the reads, so its count is now the exact input to rasusa.
        ch_read_counts = ch_dehosted_json
            .map { meta, json ->
                // The dehoster stubs emit an empty report. Treat that as "under
                // the cap" so -stub-run still exercises this path; a non-empty
                // but malformed report is a real error and throws.
                def counts = dehostCounts(json.text.trim())
                def reads_out = (counts ? counts.reads_out : 0L) as long
                // Every dehoster counts records rather than pairs, so for
                // paired-end input reads_out is twice the number of pairs.
                [meta.id, meta.single_end ? reads_out : reads_out.intdiv(2)]
            }

        ch_gated = ch_dehosted_fastq
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_read_counts)
            .branch { _id, meta, reads, count ->
                downsample: count > params.max_reads
                    return [meta, reads]
                keep: true
                    return [meta, reads]
            }

        RASUSA_READS (
            ch_gated.downsample,
            params.max_reads,
            params.downsample_seed
        )

        ch_clean_reads = ch_gated.keep.mix(RASUSA_READS.out.reads)
    }
    else {
        ch_clean_reads = ch_dehosted_fastq
    }

    // A sample under the cap bypasses rasusa and is absent here, so the report
    // falls back to the lower-priority record of wherever its reads last went.
    if (params.max_reads) {
        ch_report_records = ch_report_records.mix(
            RASUSA_READS.out.reads.map { meta, reads ->
                [ reportReads(meta.id, 'downsampling', 5, reads), [] ]
            }
        )
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

    //
    // MODULE: QQQ_RUN_REPORT -- what processed the run, what each step removed,
    // and where every read went
    //
    // One HTML page and one JSON file, published to run_report/. The JSON is the
    // contract for machines; its schema is in docs/output.md.
    //
    // The division of labour is deliberate. This workflow resolves only what
    // Nextflow alone knows -- publish paths under params.outdir, database paths
    // and whatever checksums sit beside them, the versions topic, run metadata --
    // into a manifest of plain records. Every number is parsed from the tools'
    // own reports by bin/qqq_run_report.py, where the module's nf-test covers it.
    //
    // Published paths are built from the emitted file names rather than from
    // patterns, so a module renaming its outputs cannot leave the report pointing
    // at files that do not exist. The directories are the ones conf/modules.config
    // publishes each process to; the pipeline nf-test asserts that every path in
    // the JSON exists in the test outdir, which is what keeps the two in step.
    //
    // It runs after MultiQC on purpose: it is the last thing the run produces,
    // and it records MultiQC's version, which MultiQC cannot put on the topic.
    ch_report_versions = topic_versions.versions_tuple
        .mix(MULTIQC.out.versions)
        .map { process, tool, version -> [ process: process.toString(), tool: tool.toString(), version: version.toString() ] }
        .toList()

    ch_report_inputs = ch_report_records
        .toList()
        .map { records -> [ records: records ] }
        .combine(ch_report_versions.map { versions -> [ versions: versions ] })
        .combine(MULTIQC.out.report.map { _meta, report -> [ multiqc: report.name ] })

    ch_report_manifest = ch_report_inputs
        .map { records, versions, multiqc ->
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(
                qqqReportManifest(records.records.collect { record -> record[0] }, versions.versions, multiqc.multiqc)
            ))
        }
        .collectFile(name: 'qqq_run_manifest.json')

    QQQ_RUN_REPORT (
        ch_report_manifest,
        ch_report_inputs.map { records, _versions, _multiqc -> records.records.collectMany { record -> record[1] } }
    )

    emit:
    run_report     = QQQ_RUN_REPORT.out.json                          // channel: /path/to/qqq_run_report.json
    host_bam       = ch_host_bam                                      // channel: [ val(meta), path(bam) ]
    host_bai       = ch_host_bai                                      // channel: [ val(meta), path(bai) ]
    clean_reads    = ch_clean_reads                                   // channel: [ val(meta), path(fastq) ]
    rrna_reads     = ch_rrna_reads                                    // channel: [ val(meta), path(fastq) ]
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report } // channel: /path/to/multiqc_report.html

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

//
// ---- Run report helpers -------------------------------------------------------
//
// params.outdir as the report should print it: an object-store or remote URI
// unchanged, a local path made absolute, either way without a trailing slash.
// Absolute, because a relative path in a JSON file means nothing to whatever
// reads it later from somewhere else.
//
def reportOutdir() {
    def outdir = params.outdir.toString().replaceAll('/+$', '')
    return outdir.contains('://') ? outdir : file(outdir).toAbsolutePath().toString()
}

//
// Where a process's emitted files land once published to `dir` under outdir.
// Every publishDir this pipeline uses keeps the file name (saveAs drops only
// versions.yml), so the emitted name is the published name. Takes paths, or
// plain names for run-level files that exist only after this is called.
//
def reportPublished(dir, files) {
    return [files].flatten().collect { f ->
        def name = f instanceof CharSequence ? f.toString() : f.name
        def rel = "${dir}/${name}".toString()
        [ path: "${reportOutdir()}/${rel}".toString(), relative_path: rel ]
    }
}

//
// One candidate answer to "which step last wrote this sample's reads". The
// report takes the highest priority it sees per sample; the priorities are the
// order the steps run in, so whichever step touched the reads last wins. That
// step's output is what main.nf publishes to reads/ (from ch_clean_reads), and
// nothing else publishes reads, so every candidate points there and only the
// winner's files exist. (The lowest candidate, the samplesheet input, is the
// exception: it is not a task output, so it is not published; see above.)
//
def reportReads(sample, step, priority, reads) {
    return [ kind: 'reads', sample: sample, step: step, priority: priority, published: true,
             files: reportPublished('reads', reads) ]
}

//
// The commit this run's code came from. workflow.commitId is set only when
// Nextflow pulled the pipeline itself; launched from a checkout with `nextflow
// run .`, which is how both AMDP launchers run it, it is null, so ask git.
// `dirty` is whether tracked files differ from that commit -- a run from an
// edited checkout is not reproducible from the commit alone, and should say so.
//
def reportGitCommit() {
    if (workflow.commitId) {
        return [ commit: workflow.commitId, source: 'workflow.commitId', dirty: null ]
    }
    try {
        def head = ['git', '-C', workflow.projectDir.toString(), 'rev-parse', 'HEAD'].execute()
        def sha = head.text.trim()
        head.waitFor()
        if (head.exitValue() != 0 || !sha) {
            return [ commit: null, source: null, dirty: null ]
        }
        def status = ['git', '-C', workflow.projectDir.toString(), 'status', '--porcelain', '--untracked-files=no'].execute()
        def changes = status.text.trim()
        status.waitFor()
        return [ commit: sha, source: 'git rev-parse HEAD in projectDir', dirty: status.exitValue() == 0 ? (changes as boolean) : null ]
    }
    catch (Exception _e) {
        return [ commit: null, source: null, dirty: null ]
    }
}

//
// A database or reference, described as far as is cheap. Nothing here reads a
// large file over the network: the deacon index and the rRNA reference are GB on
// S3, and hashing them at launch would cost more than the run's first tasks.
//
// So checksums come from, in order:
//   1. a sidecar beside the file: <file>.md5, .md5sum, .sha256, .sha256sum, or a
//      md5sums.txt / MD5SUMS / sha256sums.txt / SHA256SUMS in the same directory
//      listing the file by name -- a few HEAD requests on S3;
//   2. hashing it here, but only when it is on the local filesystem and no more
//      than 64 MiB, which covers test references and nothing in production.
// Otherwise they are null and checksum_note says why. Failures to inspect are
// recorded rather than thrown: a report must never fail the run it describes.
//
def reportDatabase(role, tool, location, parameters) {
    def rec = [ role: role, tool: tool, path: location.toString(), name: null, is_directory: null, size_bytes: null,
                md5: null, sha256: null, checksum_source: null, checksum_note: null, parameters: parameters ]
    try {
        def p = file(location.toString())
        rec.path = p.toUriString()
        rec.name = p.name
        rec.is_directory = p.isDirectory()
        if (!rec.is_directory) {
            rec.size_bytes = p.size()
        }
        def found = reportSidecarChecksums(p)
        if (found) {
            rec.putAll(found)
        }
        else if (rec.is_directory) {
            rec.checksum_note = 'directory; no checksum manifest found inside it'
        }
        else if (p.fileSystem == java.nio.file.FileSystems.getDefault() && rec.size_bytes <= 64 * 1024 * 1024) {
            def bytes = java.nio.file.Files.readAllBytes(p)
            rec.md5 = java.security.MessageDigest.getInstance('MD5').digest(bytes).encodeHex().toString()
            rec.sha256 = java.security.MessageDigest.getInstance('SHA-256').digest(bytes).encodeHex().toString()
            rec.checksum_source = 'computed at launch'
        }
        else {
            rec.checksum_note = 'no checksum sidecar found, and the file is remote or larger than 64 MiB so it was not hashed'
        }
    }
    catch (Exception e) {
        rec.checksum_note = "could not inspect: ${e.message}".toString()
    }
    return rec
}

def reportSidecarChecksums(p) {
    def found = [:]
    def sources = []
    def record = { hash, source ->
        def h = hash.toLowerCase()
        def algo = h ==~ /[0-9a-f]{32}/ ? 'md5' : (h ==~ /[0-9a-f]{64}/ ? 'sha256' : null)
        if (algo && !found[algo]) {
            found[algo] = h
            sources << source
        }
    }
    def candidates = []
    if (p.isDirectory()) {
        candidates = ['md5sums.txt', 'MD5SUMS', 'sha256sums.txt', 'SHA256SUMS'].collect { n -> [ p.resolve(n), null ] }
    }
    else {
        candidates = ['.md5', '.md5sum', '.sha256', '.sha256sum'].collect { ext -> [ file("${p.toUriString()}${ext}"), p.name ] }
        candidates += ['md5sums.txt', 'MD5SUMS', 'sha256sums.txt', 'SHA256SUMS'].collect { n -> [ p.parent.resolve(n), p.name ] }
    }
    candidates.each { sidecar, name ->
        try {
            if (!sidecar.exists()) {
                return
            }
            def lines = sidecar.text.readLines().findAll { line -> line.trim() }
            if (p.isDirectory()) {
                // A manifest inside a directory covers many files; point at it.
                sources << sidecar.toUriString()
                return
            }
            // `<hash>` alone, or `<hash>  [*]<name>` lines as md5sum/sha256sum write them.
            def line = lines.find { l -> l.trim().tokenize().size() == 1 }
            if (!line) {
                line = lines.find { l ->
                    def t = l.trim().tokenize()
                    t.size() > 1 && t[-1].replaceAll('^\\*', '').tokenize('/')[-1] == name
                }
            }
            if (line) {
                record.call(line.trim().tokenize()[0], sidecar.toUriString())
            }
        }
        catch (Exception _e) {
            // A sidecar we cannot read is the same as no sidecar.
        }
    }
    if (!sources) {
        return null
    }
    found.checksum_source = sources.unique().join(', ')
    return found
}

//
// The manifest QQQ_RUN_REPORT reads. Everything run-level is assembled here;
// `records` are the per-sample records mixed in as each step was wired.
//
def qqqReportManifest(records, versions, multiqc_report) {
    def git = reportGitCommit()
    def outdir = reportOutdir()

    def databases = []
    if (params.dehoster != 'none') {
        if (params.host_fasta) {
            def builder = params.dehoster == 'deacon' ? 'DEACON_INDEX' : 'BOWTIE2_BUILD (short reads); minimap2 reads the FASTA directly (long reads)'
            databases << reportDatabase('host_reference_fasta', params.dehoster, params.host_fasta,
                [ note: "index built from this FASTA at runtime by ${builder}".toString() ])
        }
        else if (params.dehoster == 'deacon') {
            databases << reportDatabase('host_index', 'deacon', params.deacon_db, [:])
        }
        else {
            databases << reportDatabase('host_index', params.dehoster, params.hostile_db, [
                index_name     : params.hostile_index,
                short_read_index: "${params.hostile_db}/${params.hostile_index}.*.bt2 (bowtie2)".toString(),
                long_read_index : "${params.hostile_db}/${params.hostile_index}.mmi (minimap2)".toString(),
            ])
        }
    }
    if (params.rrna_ref) {
        databases << reportDatabase('rrna_reference', 'bbduk', params.rrna_ref, [ k: params.rrna_k ])
    }

    // Arguments the tools do not record themselves. fastp, deacon, hostile and
    // hostess write theirs into their own reports, and the script takes them
    // from there; these come from conf/modules.config or module inputs, and are
    // labelled as such in the report because a launch-time -c override of
    // ext.args would not show here. Keep them in step with conf/modules.config.
    def declared = []
    if (params.rrna_ref) {
        declared << [ process: 'BBMAP_BBDUK', tool: 'bbduk',
                      arguments: "ref=${file(params.rrna_ref.toString()).name} k=${params.rrna_k} interleaved=f".toString() ]
    }
    if (params.max_reads) {
        declared << [ process: 'RASUSA_READS', tool: 'rasusa',
                      arguments: "--num ${params.max_reads} --seed ${params.downsample_seed}".toString() ]
    }
    if (params.dehoster == 'hostess') {
        declared << [ process: 'HOSTESS', tool: 'hostess', arguments: params.host_bam ? '--bam' : '' ]
    }

    def keep = [
        'input', 'outdir', 'publish_dir_mode', 'skip_fastp', 'save_failed_reads', 'dehoster', 'host_fasta',
        'deacon_db', 'hostile_db', 'hostile_index', 'hostile_large_input_bytes', 'host_bam', 'rrna_ref',
        'rrna_k', 'max_reads', 'downsample_seed', 'pipeline_version', 'tier',
    ]

    return [
        run: [
            pipeline        : workflow.manifest.name,
            pipeline_version: workflow.manifest.version,
            commit          : git.commit,
            commit_source   : git.source,
            commit_dirty    : git.dirty,
            revision        : workflow.revision,
            repository      : workflow.repository,
            run_name        : workflow.runName,
            session_id      : workflow.sessionId.toString(),
            started         : workflow.start.toString(),
            nextflow_version: workflow.nextflow.version.toString(),
            profile         : workflow.profile,
            container_engine: workflow.containerEngine,
            outdir          : outdir,
            launch_dir      : workflow.launchDir.toString(),
            work_dir        : workflow.workDir.toUriString(),
            project_dir     : workflow.projectDir.toString(),
            command_line    : workflow.commandLine,
            params          : keep.findAll { k -> params.containsKey(k) }.collectEntries { k ->
                def v = params[k]
                [ k, (v == null || v instanceof Boolean || v instanceof Number) ? v : v.toString() ]
            },
        ],
        settings: [
            skip_fastp       : params.skip_fastp,
            save_failed_reads: params.save_failed_reads,
            dehoster         : params.dehoster,
            host_bam         : params.host_bam,
            rrna_depletion   : params.rrna_ref as boolean,
            // The order reads pass through the removal steps; the report chains
            // each step's input to the previous step's output in this order.
            step_order       : [ 'read_filter', 'rrna_depletion', 'host_removal', 'downsampling' ],
            max_reads        : params.max_reads ?: null,
        ],
        // Run-level files, at the directories conf/modules.config publishes to.
        run_outputs: [
            run_report_json  : reportPublished('run_report', 'qqq_run_report.json')[0],
            run_report_html  : reportPublished('run_report', 'qqq_run_report.html')[0],
            multiqc_report   : reportPublished('multiqc', multiqc_report)[0],
            software_versions: reportPublished('pipeline_info', 'qqq_software_mqc_versions.yml')[0],
        ],
        databases: databases,
        versions: versions,
        declared_arguments: declared,
        records: records,
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

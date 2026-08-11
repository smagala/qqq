process HOSTESS {
    tag "${meta.id}"
    label 'process_medium'

    // hostess needs bowtie2, minimap2 and samtools, all of which the published
    // hostile image already carries, plus Python standard library only. There
    // is no hostess image yet, so this pins the internal build that adds the
    // hostess entry point to that same base.
    container 'oamd-bio-hostess:0.1.0_9890e5d_v0'

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/")
    tuple val(reference_name), path(reference_dir)

    output:
    tuple val(meta), path("*.clean*.fastq.gz")     , emit: fastq
    tuple val(meta), path("*.hostess.json")        , emit: json
    tuple val(meta), path("*.host.bam")            , emit: bam     , optional: true
    tuple val(meta), path("*.host.bam.bai")        , emit: bai     , optional: true
    tuple val(meta), path("*.host.flagstat.txt")   , emit: flagstat, optional: true
    tuple val(meta), path("*.host.stats.txt")      , emit: stats   , optional: true
    tuple val("${task.process}"), val('hostess'), eval("hostess --version | sed 's/.* //'"), topic: versions, emit: versions_hostess

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    // Sorted so mate 1 is always first: the glob that stages the reads makes
    // no ordering promise, and passing them the wrong way round would swap the
    // output mates silently.
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.simpleName }
    def reads_cmd    = meta.single_end ? "--fastq1 ${sorted_reads[0]}" : "--fastq1 ${sorted_reads[0]} --fastq2 ${sorted_reads[1]}"
    """
    hostess clean \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index ${reference_dir}/${reference_name} \\
        --prefix ${prefix} \\
        --output . \\
        --force
    """

    stub:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def want_bam = (task.ext.args ?: '').contains('--bam')
    def fastqs   = meta.single_end ? "${prefix}.clean.fastq.gz" : "${prefix}.clean_1.fastq.gz ${prefix}.clean_2.fastq.gz"
    // The BAM outputs are optional, so the stub creates them only when the
    // arguments ask for them. Otherwise a stub run would not exercise the same
    // channel shape as a real one.
    def bam_stub = want_bam
        ? "touch ${prefix}.host.bam ${prefix}.host.bam.bai ${prefix}.host.flagstat.txt ${prefix}.host.stats.txt"
        : "true"
    """
    for f in ${fastqs}; do
        echo "" | gzip > \$f
    done
    ${bam_stub}

    printf '%s\\n' \\
        '[' \\
        '    {' \\
        '        "version": "stub",' \\
        '        "aligner": "bowtie2",' \\
        '        "index": "${reference_name}",' \\
        '        "options": [],' \\
        '        "fastq1_in_name": "${prefix}",' \\
        '        "reads_in": 0,' \\
        '        "reads_out": 0,' \\
        '        "reads_removed": 0,' \\
        '        "reads_removed_proportion": 0.0' \\
        '    }' \\
        ']' > ${prefix}.hostess.json
    """
}

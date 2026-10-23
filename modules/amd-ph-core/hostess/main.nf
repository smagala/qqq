process HOSTESS {
    tag "${meta.id}"
    label 'process_medium'

    // hostess ships inside the hostile image rather than one of its own. It
    // needs bowtie2, minimap2 and samtools, which that image already carries,
    // and no third-party Python packages, so it adds well under a megabyte;
    // a separate image would only duplicate the aligners.
    //
    // No public image exists: hostess has no conda package and no BioContainer.
    // Consumers override the registry and the exact tag via their own config,
    // as with the other AMDP images.
    container 'oamd-bio-hostile:2.0.2'

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/")
    tuple val(reference_name), path(reference_dir)
    val aligner

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
    // `aligner` decides both the flag and the shape of --index, because for
    // hostess those are one decision: bowtie2 wants an index PREFIX shared by
    // several .bt2 files, minimap2 wants a single reference FASTA or .mmi FILE.
    // Taking it as an input rather than reading a meta field keeps this module
    // to meta.id and meta.single_end, and leaves the caller -- which already
    // knows whether it is routing long or short reads -- the single source of
    // truth for the choice.
    //
    // hostess's own --aligner also accepts `auto`, and this module deliberately
    // does not: choosing the shape of --index IS choosing the aligner, so `auto`
    // would mean building a bowtie2-style prefix and then letting hostess pick
    // minimap2, or the reverse. Every unrecognised value fails here rather than
    // silently taking the bowtie2 branch below.
    def index_ref = aligner == 'minimap2'
        ? (reference_dir.isDirectory() ? "${reference_dir}/${reference_name}.mmi" : reference_dir.toString())
        : "${reference_dir}/${reference_name}"
    // Sorted so mate 1 is always first: the glob that stages the reads makes
    // no ordering promise, and passing them the wrong way round would swap the
    // output mates silently.
    //
    // On the full name, not simpleName. simpleName strips every extension, so
    // `sample.R1.fastq.gz` and `sample.R2.fastq.gz` both reduce to `sample`,
    // the sort has nothing left to order on, and a pair arriving reversed stays
    // reversed -- which is the swap this line exists to prevent. Underscored
    // names like `sample_1.fastq.gz` sorted correctly under either key, which
    // is why this held up in practice.
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.name }
    // Quoted: these are staged filenames, so a space or a shell metacharacter in
    // an input name would otherwise split into separate arguments.
    def reads_cmd    = meta.single_end
        ? "--fastq1 \"${sorted_reads[0]}\""
        : "--fastq1 \"${sorted_reads[0]}\" --fastq2 \"${sorted_reads[1]}\""
    // Deliberately AFTER the declarations above, not before them. On Nextflow
    // 24.10.5 -- which this repository's nf-test matrix still covers -- a
    // statement ahead of them makes the parser read `def index_ref = aligner
    // == ...` as redeclaring an input and refuse to compile the module with
    // "Variable `aligner` already defined in the process scope". Nothing here
    // depends on the order otherwise: the declarations build strings, and this
    // throws before any of them is used.
    if (!(aligner in ['bowtie2', 'minimap2'])) {
        error("HOSTESS requires aligner to be 'bowtie2' or 'minimap2', not '${aligner}' (sample '${meta.id}'). hostess accepts 'auto', but this module has to know which index form to build.")
    }
    """
    hostess clean \\
        ${args} \\
        --aligner ${aligner} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index "${index_ref}" \\
        --prefix "${prefix}" \\
        --output . \\
        --force
    """

    stub:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def want_bam = (task.ext.args ?: '').contains('--bam')
    def fastqs   = (meta.single_end
        ? ["${prefix}.clean.fastq.gz"]
        : ["${prefix}.clean_1.fastq.gz", "${prefix}.clean_2.fastq.gz"]).collect { name -> "\"${name}\"" }.join(' ')
    // The BAM outputs are optional, so the stub creates them only when the
    // arguments ask for them. Otherwise a stub run would not exercise the same
    // channel shape as a real one.
    def bam_stub = want_bam
        ? "touch \"${prefix}.host.bam\" \"${prefix}.host.bam.bai\" \"${prefix}.host.flagstat.txt\" \"${prefix}.host.stats.txt\""
        : "true"
    // Same guard as the script path. The stub never builds an index path, so it
    // has nothing to get wrong, but a stub that accepts an aligner the real run
    // rejects is a stub that does not stand in for the real run -- and since no
    // hostess image is published, these stub tests are also the only place the
    // rejection is testable. After the declarations, for the parser reason above.
    if (!(aligner in ['bowtie2', 'minimap2'])) {
        error("HOSTESS requires aligner to be 'bowtie2' or 'minimap2', not '${aligner}' (sample '${meta.id}'). hostess accepts 'auto', but this module has to know which index form to build.")
    }
    """
    for f in ${fastqs}; do
        echo "" | gzip > "\$f"
    done
    ${bam_stub}

    printf '%s\\n' \\
        '[' \\
        '    {' \\
        '        "version": "stub",' \\
        '        "aligner": "${aligner}",' \\
        '        "index": "${reference_name}",' \\
        '        "options": [],' \\
        '        "fastq1_in_name": "${prefix}",' \\
        '        "reads_in": 0,' \\
        '        "reads_out": 0,' \\
        '        "reads_removed": 0,' \\
        '        "reads_removed_proportion": 0.0' \\
        '    }' \\
        ']' > "${prefix}.hostess.json"
    """
}

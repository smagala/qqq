process DEACON_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deacon:0.17.0--h3edb6b3_0'
        : 'quay.io/biocontainers/deacon:0.17.0--h3edb6b3_0'}"

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/")
    path index

    output:
    tuple val(meta), path("*.clean*.fq.gz"), emit: fastq
    tuple val(meta), path("*.deacon.json") , emit: json
    tuple val("${task.process}"), val('deacon'), eval("deacon --version | sed 's/.* //'"), topic: versions, emit: versions_deacon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Sorted so mate 1 is always first: the glob that stages the reads makes
    // no ordering promise, and passing them the wrong way round would swap the
    // output mates silently. Matches how hostess does it.
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.simpleName }
    def reads_cmd    = meta.single_end ? "${sorted_reads[0]}" : "${sorted_reads[0]} ${sorted_reads[1]}"
    // Named for the mate they came from, and kept distinct from the input names
    // so an input already called <prefix>.fq.gz cannot be overwritten.
    def outputs      = meta.single_end ? ["${prefix}.clean.fq.gz"] : ["${prefix}.clean_1.fq.gz", "${prefix}.clean_2.fq.gz"]
    def fifos        = outputs.collect { output -> output.replaceAll(/\.fq\.gz$/, '.fifo') }
    // -o is the first output and -O the second; deacon takes no -O at all for
    // single-end input.
    def out_flags    = ['-o', '-O']
    def gz_args      = (0..<outputs.size()).collect { i -> "${out_flags[i]} ${outputs[i]}" }.join(' ')
    def fifo_args    = (0..<fifos.size()).collect { i -> "${out_flags[i]} ${fifos[i]}" }.join(' ')
    def consumers    = (0..<outputs.size()).collect { i ->
        "bgzip --threads ${task.cpus} -c < ${fifos[i]} > ${outputs[i]} &\n        bgzip_pid_${i}=\$!"
    }.join('\n        ')
    def waits        = (0..<outputs.size()).collect { i -> "wait \$bgzip_pid_${i}" }.join('\n        ')
    """
    # Deacon writes plain gzip when handed a .gz output path. bgzip writes BGZF
    # instead: still a valid gzip stream, but block-compressed, so downstream
    # tools can index it and read it in parallel. It is also what hostess
    # produces, through samtools, so the two dehosters emit the same format.
    #
    # Not every image carries bgzip — the public BioContainer does not, the AMDP
    # image does — so it is used when present and deacon's own parallel gzip
    # writer is the fallback. The reads are identical either way.
    #
    # The bgzip path streams through FIFOs rather than writing the FASTQ out and
    # compressing it afterwards, which would put the whole uncompressed sample
    # on disk. Each consumer is waited on explicitly: a bgzip that dies would
    # otherwise leave a truncated output and a zero exit status.
    run_deacon() {
        deacon filter \\
            ${args} \\
            --threads ${task.cpus} \\
            --deplete \\
            --summary ${prefix}.deacon.json \\
            ${index} \\
            ${reads_cmd} \\
            "\$@"
    }

    if command -v bgzip > /dev/null 2>&1; then
        mkfifo ${fifos.join(' ')}
        ${consumers}
        run_deacon ${fifo_args}
        ${waits}
        rm -f ${fifos.join(' ')}
    else
        run_deacon ${gz_args}
    fi
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def outputs = meta.single_end ? ["${prefix}.clean.fq.gz"] : ["${prefix}.clean_1.fq.gz", "${prefix}.clean_2.fq.gz"]
    """
    for f in ${outputs.join(' ')}; do
        echo "" | gzip > \$f
    done

    printf '%s\\n' \\
        '{' \\
        '    "version": "deacon stub",' \\
        '    "index": "${index}",' \\
        '    "deplete": true,' \\
        '    "seqs_in": 0,' \\
        '    "seqs_out": 0,' \\
        '    "seqs_out_proportion": 0.0,' \\
        '    "seqs_removed": 0,' \\
        '    "seqs_removed_proportion": 0.0,' \\
        '    "bp_in": 0,' \\
        '    "bp_out": 0,' \\
        '    "bp_out_proportion": 0.0,' \\
        '    "bp_removed": 0,' \\
        '    "bp_removed_proportion": 0.0' \\
        '}' > ${prefix}.deacon.json
    """
}

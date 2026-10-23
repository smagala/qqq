process RASUSA_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rasusa:5.1.0--hfa8f182_0'
        : 'biocontainers/rasusa:5.1.0--hfa8f182_0'}"

    input:
    // Staged into a subdirectory so an input can never share a path with an
    // output. An input already named <prefix>.downsampled.fastq.gz would
    // otherwise be staged as the exact path this writes: rasusa creates the
    // output before it counts reads, so it truncates through the symlink and
    // destroys the upstream file, then fails on the now-empty input. Measured
    // against 5.1.0 -- a 102-byte FASTQ came back 20 bytes and zero reads.
    // Reachable by downsampling an already-downsampled sample. deacon/filter
    // and hostess stage the same way.
    tuple val(meta), path(reads, stageAs: "input_reads/*")
    val max_reads
    val seed

    output:
    tuple val(meta), path("*.downsampled.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('rasusa'), eval("rasusa --version | sed 's/.* //'"), topic: versions, emit: versions_rasusa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Sorted so mate 1 is always first. rasusa documents -o as positional
    // against the inputs -- first -o takes the first input file -- and the glob
    // that stages the reads makes no ordering promise, so an unsorted pair
    // would swap the output mates silently.
    //
    // On the full name, not simpleName: simpleName strips every extension, so
    // `sample.R1.fastq.gz` and `sample.R2.fastq.gz` both reduce to `sample` and
    // the sort becomes a no-op that leaves a reversed pair reversed. The mate
    // designator is the only thing distinguishing them and it has to stay in
    // the key.
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.name }
    def outputs = meta.single_end
        ? ["${prefix}.downsampled.fastq.gz"]
        : ["${prefix}_1.downsampled.fastq.gz", "${prefix}_2.downsampled.fastq.gz"]
    // Every path is single-quoted, with embedded quotes closed and reopened the
    // usual way. A staged name containing whitespace would otherwise word-split
    // into several arguments and rasusa would be handed paths that do not
    // exist.
    def out_args = outputs.collect { output -> "--output '" + output.toString().replace("'", "'\\''") + "'" }.join(' \\\n        ')
    def reads_cmd = sorted_reads.collect { read -> "'" + read.toString().replace("'", "'\\''") + "'" }.join(' ')
    """
    # No --strict. rasusa keeps everything and exits 0 when the input holds
    # fewer reads than --num, which is what a caller that gates on a read count
    # upstream wants: a sample that has shrunk since it was counted comes
    # through whole rather than failing the run. Pass --strict through ext.args
    # if the cap has to be met exactly.
    #
    # Compression comes from the .gz on the output names; rasusa infers the
    # format from the extension, so there is no -Z here.
    rasusa reads \\
        ${args} \\
        --num ${max_reads} \\
        --seed ${seed} \\
        ${out_args} \\
        ${reads_cmd}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = meta.single_end
        ? ["${prefix}.downsampled.fastq.gz"]
        : ["${prefix}_1.downsampled.fastq.gz", "${prefix}_2.downsampled.fastq.gz"]
    """
    for f in ${outputs.collect { output -> "'" + output.toString().replace("'", "'\\''") + "'" }.join(' ')}; do
        echo "" | gzip > "\$f"
    done
    """
}

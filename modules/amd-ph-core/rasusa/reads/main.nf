process RASUSA_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rasusa:5.1.0--hfa8f182_0'
        : 'biocontainers/rasusa:5.1.0--hfa8f182_0'}"

    input:
    tuple val(meta), path(reads)
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
    // would swap the output mates silently. Matches how deacon/filter and
    // hostess do it.
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.simpleName }
    def outputs = meta.single_end
        ? ["${prefix}.downsampled.fastq.gz"]
        : ["${prefix}_1.downsampled.fastq.gz", "${prefix}_2.downsampled.fastq.gz"]
    def out_args = outputs.collect { output -> "--output ${output}" }.join(' \\\n        ')
    def reads_cmd = sorted_reads.join(' ')
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
    for f in ${outputs.join(' ')}; do
        echo "" | gzip > \$f
    done
    """
}

process FADSO {
    tag "${meta.id}"
    label 'process_low'

    // Upstream (hkuwahara/sleaping) publishes no releases or tags, so the
    // container is built from a pinned commit and there is no conda package.
    container 'cdc-amd/fadso:20230810_1fedb59_v0'

    input:
    tuple val(meta), path(reads)
    val max_reads
    val seed

    output:
    tuple val(meta), path("*.downsampled.fastq.gz"), emit: reads
    // fadso has no --version flag and upstream publishes no releases, so the
    // pinned source commit stands in for a version per ph-core guidance. It is
    // echoed rather than hardcoded as a val so the output matches the shape
    // nf-core expects for the versions topic.
    tuple val("${task.process}"), val('fadso'), eval("echo 1fedb59"), topic: versions, emit: versions_fadso

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        fadso single \\
            ${args} \\
            -r ${seed} \\
            -k ${max_reads} \\
            -i ${reads} \\
            -o ${prefix}.downsampled.fastq.gz \\
            -z
        """
    }
    else {
        """
        fq=( ${reads} )
        fadso pair \\
            ${args} \\
            -r ${seed} \\
            -k ${max_reads} \\
            -1 "\${fq[0]}" \\
            -2 "\${fq[1]}" \\
            -a ${prefix}_1.downsampled.fastq.gz \\
            -b ${prefix}_2.downsampled.fastq.gz \\
            -z
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = meta.single_end
        ? "${prefix}.downsampled.fastq.gz"
        : "${prefix}_1.downsampled.fastq.gz ${prefix}_2.downsampled.fastq.gz"
    """
    echo "" | gzip > tmp.gz
    for out in ${outputs}; do
        cp tmp.gz "\$out"
    done
    rm tmp.gz
    """
}

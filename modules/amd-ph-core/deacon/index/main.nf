process DEACON_INDEX {
    tag "${fasta}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deacon:0.17.0--h3edb6b3_0'
        : 'quay.io/biocontainers/deacon:0.17.0--h3edb6b3_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval("deacon --version | sed 's/.* //'"), topic: versions, emit: versions_deacon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Indexing needs 2-4x the memory filtering does, which is why this is
    // process_medium where the filter step could get by with less.
    """
    deacon index build \\
        ${args} \\
        --threads ${task.cpus} \\
        --output ${prefix}.idx \\
        ${fasta}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idx
    """
}

// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HOSTILE_CLEAN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7caca3a47606de8e3460b35823193a471272aa6ab7cfafbf9aabf4615c9fa181/data'
        : 'community.wave.seqera.io/library/hostile:2.0.2--a7f5e5d341b6b94b'}"

    input:
    tuple val(meta)          , path(reads, stageAs: "input_reads/")
    tuple val(reference_name), path(reference_dir)
    val aligner

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    tuple val(meta), path('*.json')    , emit: json
    tuple val("${task.process}"), val('hostile'), eval("hostile --version"), emit: versions_hostile, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    // The aligner arrives as an input rather than being read out of meta:
    // ph-core allows only meta.id and meta.single_end, and the caller already
    // knows which aligner each alias is for -- it feeds them from separate
    // branches with their own reference channels. bowtie2 wants an index
    // prefix, minimap2 a FASTA or a .mmi, which is why the two differ here.
    def long_read = aligner == 'minimap2'
    def index_ref = long_read
    ? (reference_dir.isDirectory()
        ? "${reference_dir}/${reference_name}.mmi"
        : reference_dir.toString())
    : reference_name
    def cache_dir = long_read ? '.' : reference_dir.toString()
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.simpleName }
    def reads_cmd    = meta.single_end ? "--fastq1 ${sorted_reads[0]}" : "--fastq1 ${sorted_reads[0]} --fastq2 ${sorted_reads[1]}"
    """
    export HOSTILE_CACHE_DIR=${cache_dir}

    ## Reorder the reads for reproducibility
    ## Set offline as we never want this process to auto-download reference files as required input channel
    hostile \\
        clean \\
        ${args} \\
        --aligner ${aligner} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index ${index_ref} \\
        --output . \\
        --reorder \\
        --airplane \\
        | tee > ${prefix}.json
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def cache_dir  = aligner == 'minimap2' ? '.' : reference_dir.toString()
    def fake_read2   = !meta.single_end ? "echo '' | gzip > ${prefix}.clean_2.fastq.gz" : ""
    """
    export HOSTILE_CACHE_DIR=${cache_dir}
    echo "" | gzip > ${prefix}.clean_1.fastq.gz
    ${fake_read2}

    touch ${prefix}.json
    """
}

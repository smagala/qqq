process FALCO {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/falco:1.2.5--h077b44d_0'
        : 'biocontainers/falco:1.2.5--h077b44d_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc_report.html"), emit: html
    tuple val(meta), path("*_fastqc_data.txt")   , emit: txt
    tuple val(meta), path("*_summary.txt")       , emit: summary
    tuple val("${task.process}"), val('falco'), eval("falco --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_falco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (reads.toList().size() == 1) {
        """
        falco ${args} \\
            --threads ${task.cpus} \\
            -D ${prefix}_fastqc_data.txt \\
            -S ${prefix}_summary.txt \\
            -R ${prefix}_fastqc_report.html \\
            ${reads}
        """
    }
    else {
        // Falco writes to fixed output names, so a single invocation over both
        // mates overwrites R1's results with R2's. Run once per mate and derive
        // per-mate output names from the input filename instead.
        """
        for read_file in ${reads}; do
            read_base=\$(basename "\$read_file")
            read_base=\${read_base%.gz}
            read_base=\${read_base%.fastq}
            read_base=\${read_base%.fq}
            falco ${args} \\
                --threads ${task.cpus} \\
                -D "\${read_base}_fastqc_data.txt" \\
                -S "\${read_base}_summary.txt" \\
                -R "\${read_base}_fastqc_report.html" \\
                "\$read_file"
        done
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (reads.toList().size() == 1) {
        """
        touch ${prefix}_fastqc_data.txt
        touch ${prefix}_summary.txt
        touch ${prefix}_fastqc_report.html
        """
    }
    else {
        """
        for read_file in ${reads}; do
            read_base=\$(basename "\$read_file")
            read_base=\${read_base%.gz}
            read_base=\${read_base%.fastq}
            read_base=\${read_base%.fq}
            touch "\${read_base}_fastqc_data.txt"
            touch "\${read_base}_summary.txt"
            touch "\${read_base}_fastqc_report.html"
        done
        """
    }
}

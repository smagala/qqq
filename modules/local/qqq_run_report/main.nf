// QQQ_RUN_REPORT -- the run-level record of what qqq did, as one HTML page and one JSON file.
//
// Runs once, after everything else including MultiQC. It answers three questions a reviewer
// or a downstream script otherwise has to piece together from a dozen directories: exactly
// what processed the run (tools, versions, arguments, databases with paths and checksums),
// what each step removed from each sample, and where the reads -- kept and removed -- were
// published. The JSON schema is documented in docs/output.md; see bin/qqq_run_report.py.
//
// The workflow resolves everything only Nextflow knows (publish paths under params.outdir,
// database paths, the versions topic, run metadata) into `manifest`. The script parses the
// tools' own reports for the numbers, so that parsing is what this module's tests cover.

process QQQ_RUN_REPORT {
    tag "qqq_run_report"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // The same image viralrecon's PANEL_RUN_REPORT uses. The script needs only the Python
    // standard library, so any pinned Python 3 image would do; reusing one already in the
    // AMDP cache saves a pull. AMDP deployments point this at oamd-bio-python through
    // oamd-bio-nf-configs.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/dd746bd756e3646aa81b81a153f330e415378431b04defeade7b9ba1661aa2e2/data'
        : 'community.wave.seqera.io/library/pysam_pandas_python_pip_numpy:4420b6b7152a9e82'}"

    input:
    path manifest                              // JSON from the workflow: run, databases, versions, per-sample records
    path reports, stageAs: 'reports/*'         // tool reports the manifest names: fastp/fastplong JSON, dehoster JSON, bbduk stats

    output:
    path "*.html", emit: html
    path "*.json", emit: json
    // Not sent to the `versions` topic, for the reason MULTIQC is not: this process's input
    // waits on that topic being complete, so emitting to it would hang the pipeline. The
    // script records its own Python version in the report instead.
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'qqq_run_report'
    """
    qqq_run_report.py \\
        --manifest '${manifest}' \\
        --reports-dir reports \\
        --prefix '${prefix}' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: 'qqq_run_report'
    """
    touch '${prefix}.html'
    echo '{}' > '${prefix}.json'
    """
}

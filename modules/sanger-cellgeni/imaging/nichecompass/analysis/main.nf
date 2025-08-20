process NICHECOMPASS_ANALYSIS {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    input:
    tuple val(timestamp), path(nichecompass_dir)

    output:
    path "analysis_${prefix}_${timestamp}.ipynb", emit: notebook
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "nichecompass"  //TODO <- Is specifying default value here OK in nf-core style?
    def timestamp = task.ext.timestamp ?: ''
    """

    papermill \\
        "${moduleDir}/resources/usr/bin/nichecompass_analyse_sample_integration.ipynb" \\
        "analysis_${prefix}_${timestamp}.ipynb" \\
        -p outdir    "." \\
        -p prefix    "${prefix}" \\
        -r timestamp "${timestamp}" \\
        ${args} \\
        --kernel python3 \\
        --request-save-on-cell-execute \\
        --progress-bar \\
        --log-level INFO \\
        --log-output \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "nichecompass"
    def timestamp = task.ext.timestamp ?: ''
    """
    touch "analysis_${prefix}_${timestamp}.ipynb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

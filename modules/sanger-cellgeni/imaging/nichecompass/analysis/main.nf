process NICHECOMPASS_ANALYSIS {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    input:
    tuple val(meta), path(nichecompass_dir), path(timestamp)

    output:
    tuple val(meta), path("${nichecompass_dir}_analysis/"), emit: nichecompass_dir
    tuple val(meta), path("analysis_*.ipynb"), emit: notebook
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ts     = file(timestamp).text.trim()
    """
    rsync -a "${nichecompass_dir}" "${nichecompass_dir}_analysis"

    papermill \\
        "${moduleDir}/resources/usr/bin/nichecompass_analyse_sample_integration.ipynb" \\
        "analysis_${ts}.ipynb" \\
        -p nichecompass_dir  "${nichecompass_dir}_analysis" \\
        ${args} \\
        --kernel python3 \\
        --request-save-on-cell-execute \\
        --progress-bar \\
        --log-level INFO \\
        --log-output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rsync -a "${nichecompass_dir}" "${nichecompass_dir}_analysis"
    touch "analysis_${timestamp}.ipynb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

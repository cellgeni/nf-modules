process NICHECOMPASS_ANALYSIS {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    stageInMode 'copy'

    input:
    tuple val(meta), path(nichecompass_dir)

    output:
    tuple val(meta), path("${nichecompass_dir}"), emit: nichecompass_dir
    tuple val(meta), path("analysis_*.ipynb"), emit: notebook
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ts="\$(grep -oE '[0-9]{8}_[0-9]{6}' "${nichecompass_dir}/timestamp.txt" | head -n 1)"
    if [ -z "\$ts" ]; then
      echo "ERROR: Could not parse timestamp" >&2
      exit 1
    fi

    papermill \\
        "${moduleDir}/resources/usr/bin/nichecompass_analyse_sample_integration.ipynb" \\
        "analysis_\${ts}.ipynb" \\
        -p nichecompass_dir  "${nichecompass_dir}" \\
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
    def args = task.ext.args ?: ''
    """
    ts="\$(grep -oE '[0-9]{8}_[0-9]{6}' "${nichecompass_dir}/timestamp.txt" | head -n 1)"
    if [ -z "\$ts" ]; then
      echo "ERROR: Could not parse timestamp" >&2
      exit 1
    fi

    touch "analysis_\${ts}.ipynb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

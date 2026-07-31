process SEGTRAQ_BASELINE {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/cellgeni/segtraq:latest"

    input:
    tuple val(meta), path(zarr_dir)

    output:
    tuple val(meta), path("${prefix}_baseline_obs.csv"), emit: obs
    tuple val(meta), path("${prefix}_baseline_summary.json"), emit: summary
    tuple val("${task.process}"), val('segtraq'), eval('pip show segtraq | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segtraq_baseline.py \\
        --zarr_dir "${zarr_dir}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_baseline_obs.csv"
    echo '{"num_cells":0,"perc_unassigned_transcripts":0.0}' > "${prefix}_baseline_summary.json"
    """
}

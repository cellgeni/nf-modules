process SEGTRAQ_CLUSTERINGSTABILITY {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/segtraq:latest"
    // conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(zarr_dir)

    output:
    tuple val(meta), path("${prefix}_clustering_stability_obs.csv"),     emit: obs
    tuple val(meta), path("${prefix}_clustering_stability_metrics.json"), emit: metrics
    tuple val("${task.process}"), val('segtraq'), eval('pip show segtraq | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segtraq_clusteringstability.py \\
        --zarr_dir "${zarr_dir}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_clustering_stability_obs.csv"
    echo '{}' > "${prefix}_clustering_stability_metrics.json"
    """
}

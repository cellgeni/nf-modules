process SEGTRAQ_POINTSTATS {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/cellgeni/segtraq:latest"
    // conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(zarr_dir), path(labeled_h5ad)

    output:
    tuple val(meta), path("${prefix}_point_stats.csv"), emit: stats
    tuple val("${task.process}"), val('segtraq'), eval('pip show segtraq | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segtraq_pointstats.py \\
        --zarr_dir "${zarr_dir}" \\
        --labeled_h5ad "${labeled_h5ad}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_point_stats.csv"
    """
}

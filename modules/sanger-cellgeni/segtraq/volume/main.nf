process SEGTRAQ_VOLUME {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/cellgeni/segtraq:latest"

    input:
    tuple val(meta), path(zarr_dir)

    output:
    tuple val(meta), path("${prefix}_volume_similarity_obs.csv"), emit: obs
    tuple val("${task.process}"), val('segtraq'), eval('pip show segtraq | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segtraq_volume.py \\
        --zarr_dir "${zarr_dir}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_volume_similarity_obs.csv"
    """
}

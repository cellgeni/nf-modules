process SEGTRAQ_SUPERVISED {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/segtraq:latest"

    input:
    tuple val(meta), path(zarr_dir), path(labeled_h5ad), path(ref_h5ad)

    output:
    tuple val(meta), path("${prefix}_supervised_obs.csv"), emit: obs
    tuple val(meta), path("${prefix}_contamination_matrix.csv"), emit: contamination
    tuple val("${task.process}"), val('segtraq'), eval('pip show segtraq | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segtraq_supervised.py \\
        --zarr_dir "${zarr_dir}" \\
        --labeled_h5ad "${labeled_h5ad}" \\
        --ref_h5ad "${ref_h5ad}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_supervised_obs.csv"
    touch "${prefix}_contamination_matrix.csv"
    """
}

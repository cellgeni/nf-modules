process SEGTRAQ_LABELTRANSFER {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/segtraq:latest"
    // conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(zarr_dir), path(ref_h5ad)

    output:
    tuple val(meta), path("${prefix}_labeled.h5ad"), emit: labeled
    tuple val("${task.process}"), val('segtraq'), eval('pip show segtraq | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segtraq_labeltransfer.py \\
        --zarr_dir "${zarr_dir}" \\
        --ref_h5ad "${ref_h5ad}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_labeled.h5ad"
    """
}

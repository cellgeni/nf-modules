process RAW2FEATURES {
    tag "${meta.id}"
    label 'process_gpu'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/raw2features:0.2.0'
        : 'quay.io/cellgeni/raw2features:0.2.0'}"

    input:
    tuple val(meta), path(slide)

    output:
    tuple val(meta), path("${prefix}.embeddings.zarr"), emit: embeddings
    tuple val("${task.process}"), val('raw2features'), eval("raw2features version"), topic: versions, emit: versions_raw2features

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    raw2features embed \\
        ${slide} \\
        . \\
        ${args}

    mv *.embeddings.zarr ${prefix}.embeddings.zarr
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.embeddings.zarr/grids
    touch ${prefix}.embeddings.zarr/.zgroup
    """
}

process OMEZARR_VALIDATE {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/clij2:0.30'
        : 'quay.io/cellgeni/clij2:0.30'}"

    input:
    tuple val(meta), val(out_zarr), val(tiles_json)

    output:
    tuple val(meta), val(out_zarr), emit: out_zarr
    tuple val("${task.process}"), val('omezarr_validate'), eval("validate.py version"), topic: versions, emit: versions_omezarr_validate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    validate.py run \\
        --out_zarr ${out_zarr} \\
        --tiles_json '${tiles_json}' \\
        ${args}
    """

    stub:
    """
    echo "stub: skipping validation for ${out_zarr}"
    """
}

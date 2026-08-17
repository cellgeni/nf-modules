process OMEZARR_CREATEEMPTYCONTAINER {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/clij2:0.30'
        : 'quay.io/cellgeni/clij2:0.30'}"

    input:
    tuple val(meta), val(tiles_json), path(root_folder)

    output:
    tuple val(meta), path("${out_zarr_name}"), emit: out_zarr
    tuple val("${task.process}"), val('omezarr_createemptycontainer'), eval("create_empty_container.py version"), topic: versions, emit: versions_omezarr_createemptycontainer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_zarr_name = "${prefix}.ome.zarr"
    """
    create_empty_container.py run \\
        --root_folder ${root_folder} \\
        --tiles_json '${tiles_json}' \\
        --out_zarr ${out_zarr_name} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_zarr_name = "${prefix}.ome.zarr"
    """
    mkdir -p ${out_zarr_name}
    """
}

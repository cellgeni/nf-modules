process OMEZARR_PREPROCESS {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/clij2:0.29'
        : 'quay.io/cellgeni/clij2:0.29'}"

    input:
    tuple val(meta), path(root_folder), val(image_id), val(hcs_path), val(out_zarr)
    path psf_folder, stageAs: 'psfs'

    output:
    tuple val(meta), val(out_zarr), val(image_id), emit: fovs
    tuple val("${task.process}"), val('omezarr_preprocess'), eval("process.py version"), topic: versions, emit: versions_omezarr_preprocess

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    process.py run \\
        --root_folder ${root_folder} \\
        --out_zarr ${out_zarr} \\
        --hcs_path ${hcs_path} \\
        --psf_folder ${psf_folder} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}
    """
}

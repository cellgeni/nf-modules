process IMAGING_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/clij2:0.28"

    input:
    tuple val(meta), val(out_folder), path(root_folder), val(image_id), val(index)
    path(psf_folder, stageAs: 'psfs')

    output:
    tuple val(meta), path("${out_folder}/${out_img_name}"), emit: fovs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_img_name = "${image_id}.tif"
    """
    process.py run \\
        --root_folder ${root_folder} \\
        --index ${index} \\
        --out_img_name ${out_img_name} \\
        --out_folder ${out_folder} \\
        --psf_folder ${psf_folder} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imagepreprocess: \$(process.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_img_name = "${image_id}.tif"
    """
    mkdir -p ${out_folder}
    touch ${out_folder}/${out_img_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imagepreprocess: \$(process.py version)
    END_VERSIONS
    """
}

process SPATIAL_SCRAFTPALOM {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/scraft:latest'
        : 'quay.io/cellgeni/scraft:latest'}"

    input:
    tuple val(meta), path(ref), path(moving)

    output:
    tuple val(meta), path("${output_img_name}"), emit: aligned_moving_image
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_img_name = "aligned_" + moving.baseName.replaceFirst(/\.ome\.tif$/, '') + ".ome.tif"
    """
    scraft registration palom \\
        ${ref} ${moving} \\
        ${args} \\

    mv palom_aligned.ome.tif ${output_img_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scraft: \$(scraft info)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scraft: \$(scraft info)
    END_VERSIONS
    """
}

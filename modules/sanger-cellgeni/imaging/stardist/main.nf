process IMAGING_STARDIST {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/tiled_stardist:0.9.1-0.0.2'
        : 'quay.io/cellgeni/tiled_stardist:0.9.1-0.0.2'}"

    input:
    tuple val(meta), val(x_min), val(y_min), val(x_max), val(y_max), path(image)

    output:
    tuple val(meta), path("${output_name}"), emit: wkts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${x_min}_${y_min}_${x_max}_${y_max}"
    output_name = "${prefix}_sd_outlines.wkt"
    """
    stardist_helper.py run \\
        --image-path ${image} \\
        --x_min ${x_min} \\
        --x_max ${x_max} \\
        --y_min ${y_min} \\
        --y_max ${y_max} \\
        --output_name ${output_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(stardist_helper.py version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${x_min}_${y_min}_${x_max}_${y_max}"
    output_name = "${prefix}_sd_outlines.wkt"
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(stardist_helper.py version)
    END_VERSIONS
    """
}

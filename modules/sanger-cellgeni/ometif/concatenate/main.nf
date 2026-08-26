process OMETIF_CONCATENATE {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/scraft:latest'
        : 'quay.io/cellgeni/scraft:latest'}"

    input:
    tuple val(meta), path(images)

    output:
    tuple val(meta), path("${prefix}_hyperstack.ome.tif"), emit: hyperstack
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    concat_ome_tiffs.py \\
        ${images.join(' ')} \\
        --output ${prefix}_hyperstack.ome.tif \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concat_ome_tiffs: \$(concat_ome_tiffs.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_hyperstack.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concat_ome_tiffs: \$(concat_ome_tiffs.py --version)
    END_VERSIONS
    """
}

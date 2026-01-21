process SPATIAL_CONCATH5ADS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/scraft:latest'
        : 'quay.io/cellgeni/scraft:latest'}"

    input:
    tuple val(meta), path(anndatas)

    output:
    tuple val(meta), path("${prefix}_concatenated.h5ad"), emit: concatenated_h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    concat_h5ads.py \\
        ${anndatas.join(' ')} \\
        --output ${prefix}_concatenated.h5ad \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concat_h5ads: \$(concat_h5ads.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}_concatenated.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concat_h5ads: \$(concat_h5ads.py --version)
    END_VERSIONS
    """
}

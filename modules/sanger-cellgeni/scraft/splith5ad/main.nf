process SCRAFT_SPLITH5AD {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/scraft:latest'
        : 'quay.io/cellgeni/scraft:latest'}"

    input:
    tuple val(meta), path(merged_h5ad), val(obs_col)

    output:
    tuple val(meta), path("${prefix}/*.h5ad"), emit: h5ad_splits
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    scraft split-h5ad \\
        --input ${merged_h5ad} \\
        --output-dir ${prefix} \\
        --obs-col ${obs_col} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scraft: \$(scraft --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    mkdir ${prefix}
    touch ${prefix}/test.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scraft: \$(scraft --version)
    END_VERSIONS
    """
}

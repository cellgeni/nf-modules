process SPATIALDATA_EXPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/spatialdata-io:0.6.0'
        : 'quay.io/cellgeni/spatialdata-io:0.6.0'}"

    input:
    tuple val(meta), path(bundle)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: anndata
    tuple val(meta), path("${prefix}_mask.tif"), path("${prefix}_raw.tif"), emit: images, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_h5ad = "${prefix}.h5ad"
    """
    spatialdata_export.py \\
        ${bundle} \\
        -o ${out_h5ad} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata_export: \$(spatialdata_export --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata_export: \$(spatialdata_export --version)
    END_VERSIONS
    """
}

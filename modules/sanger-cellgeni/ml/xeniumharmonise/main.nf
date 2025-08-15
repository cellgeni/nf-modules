
process ML_XENIUMHARMONISE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch127/cellgen/cellgeni/dn10/xenium_run.sif':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(he_image), path(xenium_bundle)

    output:
    tuple val(meta), path("${output_dir}"), emit: harmonised_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_dir = "${prefix}_output"
 
    """
    harmonise_xenium.py run \\
        --he_image ${he_image} \\
        --xenium_bundle ${xenium_bundle} \\
        --output_dir ${output_dir} \\
        ${args}
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlxeniumharmonise: \$(harmonise_xenium.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_dir = "${prefix}_output"

    """
    touch ${prefix}_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlxeniumharmonise: \$(harmonise_xenium.py version)
    END_VERSIONS
    """
}

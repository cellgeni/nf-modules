process PROSEG_TO_BAYSOR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/proseg:3.2.0'
        : 'quay.io/cellgeni/proseg:3.2.0'}"

    input:
    tuple val(meta), path(sd_zarr)

    output:
    tuple val(meta), path("*-proseg-to-baysor-cell-polygons.geojson"), emit: cell_polygons
    tuple val(meta), path("*-proseg-to-baysor-transcript-metadata.csv"), emit: transcript_metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proseg-to-baysor \\
        ${sd_zarr} \\
        --output-transcript-metadata ${prefix}-proseg-to-baysor-transcript-metadata.csv \\
        --output-cell-polygons ${prefix}-proseg-to-baysor-cell-polygons.geojson \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-proseg-to-baysor-transcript-metadata.csv
    touch ${prefix}-proseg-to-baysor-cell-polygons.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: "stub"
    END_VERSIONS
    """
}

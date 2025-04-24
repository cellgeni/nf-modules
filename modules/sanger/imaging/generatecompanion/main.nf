process IMAGINGGENERATECOMPANION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/imagetileprocessor:0.1.13"

    input:
    tuple val(meta), path(file_with_ome_md)

    output:
    tuple val(meta), path("${out_name}"), emit: csv
    tuple val(meta), path("${out_companion}"), emit: companion
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_image_tiles.csv"
    out_companion = "${prefix}_image_companion.xml"
    """
    generate_companion.py run \\
        --file_with_ome_md ${file_with_ome_md} \\
        --tiles_csv ${out_name} \\
        --companion_xml ${out_companion} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imagepreprocess: \$(generate_companion.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_image_tiles.csv"
    out_companion = "${prefix}_image_companion.xml"
    """
    touch ${out_name}
    touch ${out_companion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imagepreprocess: \$(generate_companion.py version)
    END_VERSIONS
    """
}

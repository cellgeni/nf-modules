process IMAGING_GENERATECOMPANION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/imagetileprocessor:0.1.13"

    input:
    tuple val(meta), val(round), path(image_root_folder)

    output:
    tuple val(meta), val(out_folder), path("${out_folder}/${out_name}"), emit: csv
    tuple val(meta), val(out_folder), path("${out_folder}/${out_companion}"), emit: companion
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_image_tiles.csv"
    out_companion = "${prefix}_image.companion.ome"
    out_folder = "image${round}"
    """
    generate_companion.py run \\
        --image_root_folder ${image_root_folder} \\
        --tiles_csv ${out_name} \\
        --companion_xml ${out_companion} \\
        --out_folder ${out_folder} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imagepreprocess: \$(generate_companion.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_image_tiles.csv"
    out_companion = "${prefix}_image.companion.ome"
    out_folder = "image${round}"
    """
    mkdir -p ${out_folder}
    touch ${out_folder}/${out_name}
    touch ${out_folder}/${out_companion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imagepreprocess: \$(generate_companion.py version)
    END_VERSIONS
    """
}

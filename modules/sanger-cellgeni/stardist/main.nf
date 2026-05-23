process STARDIST {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/stardist:0.9.2-tileprocessor-0.2.1'
        : 'quay.io/cellgeni/stardist:0.9.2-tileprocessor-0.2.1'}"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("${output_name}"), emit: geojson
    tuple val("${task.process}"), val('stardist'), eval("stardist_full_image_helper.py version"), topic: versions, emit: versions_stardist

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_sd_outlines.geojson"
    """
    stardist_full_image_helper.py download-model \\
        ${args}

    stardist_full_image_helper.py run \\
        --image-path "${image}" \\
        --output-name "${output_name}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_sd_outlines.geojson"
    """
    printf '%s\\n' '{"type":"FeatureCollection","features":[]}' > "${output_name}"
    """
}

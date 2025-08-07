process IMAGING_GENERATECOMPANIONFROMFILES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/imagetileprocessor:0.1.16'
        : 'quay.io/cellgeni/imagetileprocessor:0.1.16'}"

    input:
    tuple val(meta), path(images, stageAs: "images/*"), val(file_regex)

    output:
    tuple val(meta), path("${out_companion_name}"), emit: companion
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_companion_name = "${prefix}.companion.ome"
    """
    generate_companion.py run \\
        -input_folder images \\
        -regex '${file_regex}' \\
        -out_companion_xml "${out_companion_name}" \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(generate_companion.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_companion_name = "${prefix}.companion.ome"
    """
    echo ${args}
    
    touch ${out_companion_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(generate_companion.py version)
    END_VERSIONS
    """
}

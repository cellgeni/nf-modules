process IMAGING_CONCATENATEWKTS {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/bioinfotongli/tiled_spotiflow:0.5.2'
        : 'quay.io/bioinfotongli/tiled_spotiflow:0.5.2'}"

    input:
    tuple val(meta), path(wkts)

    output:
    tuple val(meta), path("${output_name}"), emit: concatenated_peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_merged_peaks.csv"
    """
    merge_wkts.py run \\
        -output_name ${output_name} \\
        ${wkts} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(merge_wkts.py version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_merged_peaks.csv"
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(merge_wkts.py version)
    END_VERSIONS
    """
}

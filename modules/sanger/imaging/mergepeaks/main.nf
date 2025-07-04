process IMAGING_MERGEPEAKS {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/imagetileprocessor:0.1.13'
        : 'quay.io/cellgeni/imagetileprocessor:0.1.13'}"

    input:
    tuple val(meta), val(ch_ind), path(csvs)

    output:
    tuple val(meta), path("${output_name}"), emit: merged_peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_merged_peaks_ch_${ch_ind}.wkt"
    """
    merge-peaks run \
        --output_name ${output_name} \
        ${csvs} \
        ${args} \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(merge-peaks version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_merged_peaks_ch_${ch_ind}.wkt"
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(merge-peaks version)
    END_VERSIONS
    """
}

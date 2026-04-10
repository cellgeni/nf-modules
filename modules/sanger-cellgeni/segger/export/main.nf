process SEGGER_EXPORT {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/segger:v2"

    input:
    tuple val(meta), path(segmentation), path(source_dir)

    output:
    tuple val(meta), path("${out_dir}"), emit: export_dir
    tuple val("${task.process}"), val('segger'), eval('pip show segger | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger_export"
    """
    segger export \\
        -s ${segmentation} \\
        -i ${source_dir} \\
        -o ${out_dir} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger_export"
    """
    mkdir -p ${out_dir}
    """
}

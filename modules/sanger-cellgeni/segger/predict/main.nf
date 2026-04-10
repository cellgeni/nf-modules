process SEGGER_PREDICT {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/cellgeni/segger:v2"

    input:
    tuple val(meta), path(checkpoint), path(input_dir)

    output:
    tuple val(meta), path("${out_dir}/segger_segmentation.parquet"), emit: segmentation
    tuple val("${task.process}"), val('segger'), eval('pip show segger | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger"
    """
    segger predict \\
        -c ${checkpoint} \\
        -i ${input_dir} \\
        -o ${out_dir} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger"
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/segger_segmentation.parquet
    """
}

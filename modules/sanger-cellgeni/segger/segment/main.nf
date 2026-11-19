process SEGGER_SEGMENT {
    tag "${meta.id}"
    label 'process_gpu'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/segger:0.3.1'
        : 'quay.io/cellgeni/segger:0.3.1'}"

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("${out_dir}/segger_segmentation.parquet"), emit: segmentation
    tuple val("${task.process}"), val('segger'), eval('python -c "from importlib.metadata import version; print(version(\'segger\'))"'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger"
    """
    segger segment \\
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

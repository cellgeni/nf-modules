process LOCALMIRROR {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(input_data)

    output:
    tuple val(meta), path("mirrored/${mirror_name}"), emit: mirrored
    tuple val("${task.process}"), val('coreutils'), eval("cp --version | head -n 1 | sed 's/^cp (GNU coreutils) //'"), topic: versions, emit: versions_localmirror

    when:
    task.ext.when == null || task.ext.when

    script:
    mirror_name = task.ext.prefix ?: input_data.getName()
    """
    mkdir -p mirrored
    cp -rL "${input_data}" "mirrored/${mirror_name}"
    """

    stub:
    mirror_name = task.ext.prefix ?: input_data.getName()
    """
    mkdir -p mirrored
    touch "mirrored/${mirror_name}"
    """
}

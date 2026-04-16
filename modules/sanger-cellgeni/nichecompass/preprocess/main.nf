process NICHECOMPASS_PREPROCESS {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/nichecompass:0.3.1"

    input:
    tuple val(meta), path(h5ad, stageAs: "inputs/*")

    output:
    tuple val(meta), path("${prefix}_preprocessed.h5ad"), path("data"), emit: preprocessed
    tuple val("${task.process}"), val('nichecompass'), eval('pip show nichecompass | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=/tmp
    export XDG_CONFIG_HOME=/tmp/.config
    nichecompass_preprocess.py \\
        --batches ${h5ad} \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_preprocessed.h5ad"
    mkdir -p data
    """
}

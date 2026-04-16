process NICHECOMPASS_TRAIN {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/cellgeni/nichecompass:0.3.1"

    input:
    tuple val(meta), path(preprocessed_h5ad), path(nichecompass_data)

    output:
    tuple val(meta), path("${out_name}/artifacts/model"), path("${out_name}/run_config.json"), path("${out_name}/data"), emit: nichecompass_model
    tuple val(meta), path("${out_name}/artifacts/figures"), emit: nichecompass_figures, optional: true
    tuple val("${task.process}"), val('nichecompass'), eval('pip show nichecompass | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_nichecompass_dir"
    """
    export HOME=/tmp
    export XDG_CONFIG_HOME=/tmp/.config
    nichecompass_train.py \\
        --preprocessed_h5ad "${preprocessed_h5ad}" \\
        --model_dir "${out_name}" \\
        ${args}
    cp -r "${nichecompass_data}" "${out_name}/data"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_nichecompass_dir"
    """
    mkdir -p "${out_name}/artifacts/model"
    mkdir -p "${out_name}/artifacts/figures"
    mkdir -p "${out_name}/data"
    touch "${out_name}/artifacts/model/trained.h5ad"
    touch "${out_name}/run_config.json"
    """
}

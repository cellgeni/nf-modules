process NICHECOMPASS_TRAINING {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    input:
    tuple val(meta), path(h5ad, stageAs: "inputs/*")

    output:
    tuple val(meta), path("${prefix}_*"), emit: nichecompass_dir  // TODO is parsing path with glob fine?
    tuple val(meta), path("${prefix}_*/timestamp.txt"), emit: timestamp
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nichecompass_train_sample_integration.py \\
        --batches ${h5ad} \\
        --prefix "${prefix}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def timestamp = new Date().format("yyyyMMdd_HHmmss")
    """
    mkdir -p "${prefix}_${timestamp}"
    mkdir -p "${prefix}_${timestamp}/artifacts/sample_integration/nichecompass/figures"
    mkdir -p "${prefix}_${timestamp}/artifacts/sample_integration/nichecompass/models"
    mkdir -p "${prefix}_${timestamp}/data"
    touch "${prefix}_${timestamp}/artifacts/sample_integration/nichecompass/models/model.h5ad"
    touch "${prefix}_${timestamp}/train.log"
    touch "${prefix}_${timestamp}/timestamp.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

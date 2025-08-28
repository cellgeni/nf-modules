process SPATIAL_NICHECOMPASSTRAINING {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    input:
    tuple val(meta), path(h5ad, stageAs: "inputs/*")

    output:
    tuple val(meta), path("nichecompass_*"), emit: nichecompass_model
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nichecompass_train_sample_integration.py \\
        --batches ${h5ad} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    def timestamp = new Date().format("yyyyMMdd_HHmmss")
    """
    mkdir -p "nichecompass_${timestamp}"
    mkdir -p "nichecompass_${timestamp}/artifacts/figures"
    mkdir -p "nichecompass_${timestamp}/artifacts/models"
    mkdir -p "nichecompass_${timestamp}/data"
    touch "nichecompass_${timestamp}/artifacts/models/model.h5ad"
    touch "nichecompass_${timestamp}/train.log"
    touch "nichecompass_${timestamp}/timestamp.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

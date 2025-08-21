process NICHECOMPASS_TRAINING {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    input:
    tuple val(meta), path(h5ad, stageAs: "inputs/*")

    output:
    path "nichecompass_*", emit: nichecompass_outdir  // TODO is parsing path with glob fine?
    path "nichecompass_*/timestamp.txt", emit: timestamp
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nichecompass_train_sample_integration.py \\
        --batches "${h5ad}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def timestamp = new Date().format("yyyyMMdd_HHmmss")
    """
    mkdir -p "nichecompass_${prefix}_${timestamp}"
    mkdir -p "nichecompass_${prefix}_${timestamp}/artifacts/sample_integration/nichecompass/figures"
    mkdir -p "nichecompass_${prefix}_${timestamp}/artifacts/sample_integration/nichecompass/models"
    mkdir -p "nichecompass_${prefix}_${timestamp}/data"
    touch "nichecompass_${prefix}_${timestamp}/artifacts/sample_integration/nichecompass/models/model.h5ad"
    touch "nichecompass_${prefix}_${timestamp}/train.log"
    echo "${timestamp}" > "nichecompass_${prefix}_${timestamp}/timestamp.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

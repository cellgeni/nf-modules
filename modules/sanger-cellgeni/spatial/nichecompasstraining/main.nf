process SPATIAL_NICHECOMPASSTRAINING {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    input:
    tuple val(meta), path(h5ad, stageAs: "inputs/*"), val(suffix)

    output:
    tuple val(meta), path("${out_name}/artifacts/model"), path("${out_name}/run_config.json"), path("${out_name}/data"), emit: nichecompass_model
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_nichecompass_dir"
    """
    nichecompass_train_sample_integration.py \\
        --batches ${h5ad} \\
        --model_dir "${out_name}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_nichecompass_dir"
    """
    mkdir -p "${out_name}"
    mkdir -p "${out_name}/artifacts/model"
    mkdir -p "${out_name}/data"
    touch "${out_name}/artifacts/model/model.h5ad"
    touch "${out_name}/run_config.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
    END_VERSIONS
    """
}

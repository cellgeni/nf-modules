process DOTPY {
    tag "${meta.id}"
    label params.dotpy_use_gpu ? 'process_gpu' : 'process_normal'

    container "quay.io/cellgeni/dotpy:latest"

    input:
    tuple val(meta), path(sc_adata), path(sp_adata)

    output:
    tuple val(meta), path("*_weights.csv"), emit: weights
    tuple val(meta), path("*_annotations.csv"), emit: annotations
    tuple val("${task.process}"), val('dotpy'), eval('python -c "import dotpy; print(dotpy.__version__)"'), topic: versions, emit: versions_dotpy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_dot_cli.py \\
        --ref ${sc_adata} \\
        --spatial ${sp_adata} \\
        --output ${prefix} \\
        --output-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_weights.csv
    touch ${prefix}_annotations.csv
    """
}

process TACCO_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/tacco:0.5.0'
        : 'quay.io/cellgeni/tacco:0.5.0'}"

    input:
    // spatial: AnnData to annotate; reference: annotated single-cell AnnData. Both with raw counts.
    tuple val(meta), path(spatial, stageAs: 'input/*'), path(reference, stageAs: 'reference/*')

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}_annotation.csv.gz"), emit: annotation
    tuple val("${task.process}"), val('tacco'), eval('python -c "from importlib.metadata import version; print(version(\'tacco\'))"'), emit: versions_tacco, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // The reference .obs column to transfer (--annotation-key), the method (--method OT|nnls|
    // projection|svm|..., default OT), counts layers (--spatial-layer, --reference-layer) and the
    // tacco.tl.annotate parameters (--bisections, --multi-center, --max-annotation, --n-hvg,
    // --method-kwargs '{"epsilon": 0.01}', ...) go in ext.args. See tacco_annotate.py --help.
    """
    export MPLCONFIGDIR=\$PWD/.cache/matplotlib
    export NUMBA_CACHE_DIR=\$PWD/.cache/numba
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export NUMBA_NUM_THREADS=${task.cpus}

    tacco_annotate.py \\
        --spatial ${spatial} \\
        --reference ${reference} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.h5ad
    echo "" | gzip > ${prefix}_annotation.csv.gz
    """
}

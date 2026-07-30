process TANGRAM {
    tag "${meta.id}"
    label params.tangram_use_gpu ? 'process_gpu' : 'process_normal'

    container "quay.io/cellgeni/tangram:1.0.4"

    input:
    tuple val(meta), path(sc_data), path(sp_data), path(markers)

    output:
    tuple val(meta), path("tangram_aligned.h5ad"), emit: mapped
    tuple val(meta), path("tangram_gene_proj.h5ad"), emit: geneproj
    tuple val(meta), path("spatial_with_tangram_celltypes.h5ad"), emit: spatial_ct
    path "figures/*.png", optional: true, emit: figures
    path "tangram.log", emit: log
    tuple val("${task.process}"), val('tangram'), eval('python3 -c "import tangram; print(tangram.__version__)"'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def markers_arg = markers.name != 'NO_FILE' && markers.size() > 0 ? "--markers ${markers}" : ""
    """
    mkdir -p figures

    tangram_wrapper.py \\
        --scrna ${sc_data} \\
        --spatial ${sp_data} \\
        --outdir . \\
        ${markers_arg} \\
        --save-props-csv \\
        ${args}
    """

    stub:
    """
    touch tangram_aligned.h5ad tangram_gene_proj.h5ad spatial_with_tangram_celltypes.h5ad tangram.log
    mkdir -p figures
    """
}

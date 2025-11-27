process ADATA_ATTACHCELLMETADATA {
    tag "Attaching cell metadata to .obs section of AnnData object for ${meta.id}"
    container "${ workflow.containerEngine == 'singularity' ? 'docker://quay.io/cellgeni/metacells-python:latest': 'quay.io/cellgeni/metacells-python:latest' }"
    
    input:
    tuple val(meta), path(h5ad, name: 'input/*')
    path(metadata)
    
    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: '--barcode_column obs_names'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    attach_annotation.py \
        ${args} \
        --h5ad_file ${h5ad} \
        --sample_id ${meta.id} \
        --metadata ${metadata} \
        --output ${prefix}_meta.h5ad
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}

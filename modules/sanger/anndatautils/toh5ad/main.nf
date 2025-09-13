

process ANNDATAUTILS_TOH5AD {
    tag "Converting ${sample_id}'s file to .h5ad"
    container 'docker://quay.io/cellgeni/metacells-python:latest'
    
    input:
        tuple val(sample_id), path(input, name: 'input/*')
        val(delimiter)
    output:
        tuple val(sample_id), path("${sample_id}.h5ad"), emit: h5ad
        path "versions.yml", emit: versions
    script:
        """
        convert_to_h5ad.py \
            --input "${input}" \
            --sample_id "${sample_id}" \
            --delimiter "${delimiter}" \
            --output "${sample_id}.h5ad"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            anndata: \$( python -c "import anndata; print(anndata.__version__)" )
            scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
        END_VERSIONS
        """
    stub:
        """
        touch "${sample_id}.h5ad"
        touch versions.yml
        """
}

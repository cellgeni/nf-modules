process SPATIAL_TERRA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // ? '/nfs/cellgeni/singularity/images/terra.sif'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? '/nfs/cellgeni/singularity/images/terra.sif'
        : '/nfs/cellgeni/singularity/images/terra.sif'}"

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path("${out_folder}/${baseName}_processed.h5ad"), emit: processed_anndata
    tuple val(meta), path("${out_folder}/${baseName}_tokenized_data/*"), emit: tokenized_data, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    baseName = anndata.baseName
    out_folder = "${baseName}"
    """
    export NUMBA_CACHE_DIR="/tmp/numba_cache"
    run_terra.py \\
        --h5ad_path ${anndata} \\
        --tokenized_data_path ${out_folder}/${baseName}_tokenized_data \\
        --output_file ${out_folder}/${baseName}_processed.h5ad \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatial: \$(run_terra.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    out_folder = "${anndata.baseName}"
    """
    echo ${args}
    export NUMBA_CACHE_DIR="/tmp/numba_cache"
    
    mkdir -p ${out_folder}/${baseName}_tokenized_data
    touch ${out_folder}/${baseName}_processed.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatial: \$(run_terra.py --version)
    END_VERSIONS
    """
}

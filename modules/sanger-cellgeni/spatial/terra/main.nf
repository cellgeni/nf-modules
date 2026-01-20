process SPATIAL_TERRA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // ? '/nfs/cellgeni/singularity/images/terra.sif'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? '/home/ubuntu/deleteme/terra.sif'
        : '/nfs/cellgeni/singularity/images/terra.sif'}"

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path("${out_folder}/*.h5ad"), emit: processed_anndata
    tuple val(meta), path("${out_folder}/tokenized_data/*"), emit: tokenized_data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_folder = "${prefix}"
    """
    export NUMBA_CACHE_DIR="/tmp/numba_cache"
    run_terra.py \\
        --h5ad_path ${anndata} \
        --tokenized_data_path ${out_folder}/tokenized_data \
        --output_file ${out_folder}/spatial_data_processed.h5
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatial: \$(run_terra.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_folder = "${prefix}"
    """
    echo ${args}
    export NUMBA_CACHE_DIR="/tmp/numba_cache"
    
    mkdir -p ${out_folder}/tokenized_data
    touch ${out_folder}/spatial_data_processed.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatial: \$(run_terra.py --version)
    END_VERSIONS
    """
}


process ML_XENIUMHARMONISE {
    tag "$meta"

    container '/nfs/cellgeni/singularity/images/palom_tiledb.sif'

    input:
    tuple val(meta), path(he_image), path(xenium_bundle)

    output:
    tuple val(meta), path(output_dir), emit: harmonised_output
    path "versions.yml"               , emit: versions

    script:

    output_dir = "${meta}_harmonised"
 
    """
    harmonise_xenium.py main\\
        --hematoxylin_eosin_image_uri ${he_image} \\
        --xenium_bundle_uri ${xenium_bundle} \\
        --harmonised_dataset_uri ${output_dir} \\

        
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlxeniumharmonise: 1.0.0
    END_VERSIONS
    """

    stub:
    output_dir = "${meta}_harmonised"

    """
    touch ${output_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlxeniumharmonise: 1.0.0
    END_VERSIONS
    """
}


process ML_XENIUMHARMONISE {
    tag "$meta"

    container '/nfs/cellgeni/singularity/images/palom_tiledb.sif'
        

    input:
    tuple val(meta), path(he_image), path(xenium_bundle)

    output:
    tuple val(meta), path(out_dir), emit: harmonised_output
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    out_dir = "${meta}_harmonised"
 
    """
        harmonise_xenium.py main\\
            --hematoxylin_eosin_image_uri ${he_image} \\
            --xenium_bundle_uri ${xenium_bundle} \\
            --harmonised_dataset_uri ${out_dir} 
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mlxeniumharmonise: 1.0.0
        END_VERSIONS
    """

    stub:
    out_dir = "${meta}_harmonised"

    """
    touch ${out_dir}
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mlxeniumharmonise: 1.0.0
        END_VERSIONS
    """
}

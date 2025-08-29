
process ML_XENIUMHARMONISE {
    tag "$meta"

    container '/nfs/cellgeni/singularity/images/palom_tiledb.sif'
        

    input:
    tuple val(meta), path(he_image), path(xenium_bundle)

    output:
    tuple val(meta), path(out_dir), emit: harmonised_output
    path "versions.yml"                        , emit: versions

    script:
    out_dir = "${meta}_harmonised"
 
    """
        harmonise_xenium.py main\\
            --hematoxylin_eosin_image_uri ${he_image} \\
            --xenium_bundle_uri ${xenium_bundle} \\
            --harmonised_dataset_uri ${out_dir} 

        echo "v1.0.0" > versions.yml
    """

    stub:
    out_dir = "${meta}_harmonised"

    """
    touch ${out_dir}
    touch versions.yml
    """
}

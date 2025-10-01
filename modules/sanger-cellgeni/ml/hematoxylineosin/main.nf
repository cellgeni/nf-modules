process ML_HEMATOXYLINEOSIN {
    tag "$meta"
    label 'process_medium'

    input:
    tuple val(meta), path(harmonised_output)

    output:
    tuple val(meta), path(output_dir), emit: data_array
    path "versions.yml"              , emit: versions


    script:
    output_dir = "${meta}_hematoxylin_eosin.zarr"
 
    """
        
        feature_hematoxylin_eosin_image_tiles.py generate_data_array \\
            --harmonised_dataset_absolute_path ${harmonised_output} \\
            --resolution_um_per_px 0.5 \\
            --data_array_absolute_path ${output_dir} \\
        
        echo "v1.0.0" > versions.yml
    """

    stub:
    output_dir = "${meta}_hematoxylin_eosin.zarr"
    
    """
    touch ${output_dir}
    touch versions.yml
    """
}
 
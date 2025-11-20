process ML_SEGMENTATIONCELLVECTORFIELDTILES {
    tag "$meta"
    label 'process_medium'

    
    input:
    tuple val(meta), path(harmonisexenium_output)

    output:
   
    tuple val(meta), path(output_dir), emit: data_array
    path "versions.yml"              , emit: versions

    script:
    output_dir = "${meta}_segmentation_cell_vector_field_tiles.zarr"

    """
        feature_segmentation_cell_vector_field_tiles.py generate_data_array \\
            --harmonised_dataset_absolute_path ${harmonisexenium_output} \\
            --data_array_absolute_path ${output_dir} \\
        
        echo "v1.0.0" > versions.yml
    """

    stub:
    output_dir = "${meta}_segmentation_cell_vector_field_tiles.zarr"
    """
    touch ${output_dir}
    touch versions.yml
    """
}

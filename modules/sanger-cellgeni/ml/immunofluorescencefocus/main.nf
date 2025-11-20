process ML_IMMUNOFLUORESCENCEFOCUS {
    tag "$meta"
    label 'process_medium'

    input:
    tuple val(meta), path(harmonisexenium_output)

    output:
    tuple val(meta), path(out_dir), emit: data_array
    path "versions.yml"           , emit: versions


    script:
    out_dir ="${meta}_immunofluorescence_focus"
 
    """
        feature_morphology_immunofluorescence_image_tiles.py generate_data_array \\
            --harmonised_dataset_absolute_path ${harmonisexenium_output} \\
            --data_array_absolute_path ${out_dir} \\

        echo "v1.0.0" > versions.yml
    """

    stub:

    out_dir ="${meta}_immunofluorescence_focus"

    """
    touch ${out_dir}
    touch versions.yml

    """
}



process ANNDATAUTILS_TOH5AD {
    tag "Converting ${sample_id}'s file to .h5ad"
    input:
        tuple val(sample_id), path(input, name: 'input/*')
        val(delimiter)
    output:
        tuple val(sample_id), path("${sample_id}.h5ad")
    script:
        """
        convert_to_h5ad.py \
            --input ${input} \
            --sample_id ${sample_id} \
            --delimiter ${delimiter} \
            --output ${sample_id}.h5ad
        """
}

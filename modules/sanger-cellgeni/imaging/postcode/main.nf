process IMAGING_POSTCODE {
    tag "${meta.id}"

    container "quay.io/cellgeni/postcode:0.2.0"

    input:
    tuple val(meta), file(profile), file(starfish_codebook), file(barcode_0123), file(spot_loc)

    output:
    tuple val(meta), path("${out_name}"), emit: decoded_peaks
    tuple val(meta), path("${prefix}_decode_out_parameters.pickle"), optional: true
    path "versions.yml", emit: versions

    script:
    prefix = meta.id ?: "none"
    out_name = "${prefix}_decoded_spots.csv"
    def args = task.ext.args ?: ""
    """
    decode.py run \
        --spot_profile_p ${profile} \
        --spot_locations_p ${spot_loc} \
        --starfish_codebook_p ${starfish_codebook} \
        --barcode_0123_p ${barcode_0123} \
        --out_name ${out_name} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(decode.py version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}_decoded_spots.csv"
    """
    touch ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(decode.py version)
    END_VERSIONS
    """
}

process IMAGING_POSTCODEPREP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/cellgeni/postcode:0.2.0"

    input:
    tuple val(meta), path(profile), path(tabular_codebook), path(readout_file), val(R)

    output:
    tuple val(meta), path(reformatted_profile), path(starfish_codebook), path(barcode_0123_str), emit: for_decoding
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readout = readout_file ? "--readouts ${readout_file}" : ""
    reformatted_profile = "${prefix}_reformatted_profile.npy"
    starfish_codebook = "${prefix}_starfish_codebook.json"
    barcode_0123_str = "${prefix}_barcodes_0123_str.txt"
    """
    postcodeprep.py \\
        --profile ${profile} \\
        --tabular_codebook ${reformatted_profile} \\
        --out_starfish_codebook ${starfish_codebook} \\
        --out_reformatted_profile ${reformatted_profile} \\
        --out_barcodes_0123_str ${barcode_0123_str} \\
        --R ${R} \\
        ${readout} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postcodeprep: \$(postcodeprep.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    reformatted_profile = "${prefix}_reformatted_profile.npy"
    starfish_codebook = "${prefix}_starfish_codebook.json"
    barcode_0123_str = "${prefix}_barcodes_0123_str.txt"
    """
    echo ${args}
    
    touch ${reformatted_profile}
    touch ${starfish_codebook}
    touch ${barcode_0123_str}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postcodeprep: \$(postcodeprep.py --version)
    END_VERSIONS
    """
}

process PE2OMETIF {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/imagetileprocessor:0.2.0'
        : 'quay.io/cellgeni/imagetileprocessor:0.2.0'}"

    input:
    tuple val(meta), path(image_dir), val(index_name)

    output:
    tuple val(meta), path("${prefix}/${prefix}*.ome.tif"), path("${prefix}/${prefix}.companion.ome"), emit: ome_tif
    tuple val(meta), path("${prefix}/${prefix}_ch*_ffp.tiff"), optional: true, emit: ffp_maps
    tuple val(meta), path("${prefix}/${prefix}_ch*_dfp.tiff"), optional: true, emit: dfp_maps
    tuple val("${task.process}"), val('pe2ometif'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_pe2ometif

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def round_index = meta.round_index ? "_${meta.round_index}" : "_1"
    prefix = task.ext.prefix ?: "${meta.id}${round_index}"
    """
    pe_to_ome_tif.py \\
        ${image_dir}/${index_name} \\
        --output-dir ${prefix} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def round_index = meta.round_index ? "_${meta.round_index}" : "_1"
    prefix = task.ext.prefix ?: "${meta.id}${round_index}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/A01_F001_maxproj.ome.tif
    touch ${prefix}/${prefix}.companion.ome
    touch ${prefix}/${prefix}_ch01_ffp.tiff
    touch ${prefix}/${prefix}_ch01_dfp.tiff
    """
}

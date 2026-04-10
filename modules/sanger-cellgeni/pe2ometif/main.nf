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
    tuple val(meta), path("${prefix}/${prefix}*.ome.tif"), emit: ome_tif
    tuple val(meta), path("${prefix}/${prefix}.companion.ome"), optional: true, emit: companion
    tuple val("${task.process}"), val('pe2ometif'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_pe2ometif

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pe_to_ome_tif.py \\
        ${image_dir}/${index_name} \\
        --output-dir ${prefix} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/A01_F001_maxproj.ome.tif

    if [[ ! " ${args} " =~ " --no-companion " ]]; then
        _img_prefix=\$(echo "${args}" | sed -n 's/.*--prefix[[:space:]]\\+\\([^[:space:]]\\+\\).*/\\1/p')
        _img_prefix=\${_img_prefix:-plate}
        touch ${prefix}/\${_img_prefix}.companion.ome
    fi
    """
}

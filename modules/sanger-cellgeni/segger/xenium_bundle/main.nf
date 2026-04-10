process SEGGER_XENIUM_BUNDLE {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/cellgeni/segger:v2"

    input:
    tuple val(meta), path(xenium_dir), path(export_dir)

    output:
    tuple val(meta), path("${out_dir}"), emit: xenium_bundle
    tuple val("${task.process}"), val('segger'), eval('pip show segger | grep "^Version" | sed \'s/Version: //\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_xenium_segger"
    """
    # Duplicate the xenium bundle locally, resolving any symlinks from Nextflow staging
    cp -rL ${xenium_dir} ${out_dir}

    # Inject all segger export files into the duplicated bundle
    cp -r ${export_dir}/* ${out_dir}/
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_xenium_segger"
    """
    mkdir -p ${out_dir}
    """
}

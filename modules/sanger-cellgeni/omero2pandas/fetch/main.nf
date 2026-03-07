process OMERO2PANDAS_FETCH {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/omero2pandas:latest'
        : 'quay.io/cellgeni/omero2pandas:latest'}"

    secret 'OMERO_USER'
    secret 'OMERO_PASSWORD'

    input:
    tuple val(meta), val(annotation_id), val(linking_index)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}/annotation_*.csv"), emit: annotations
    tuple val("${task.process}"), val('omero2pandas_fetch'), eval("omero_download_annotations.py --version"), topic: versions, emit: versions_omero2pandas_fetch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_dir = "${prefix}"
    """
    omero_download_annotations.py \\
        --user \$OMERO_USER \\
        --password \$OMERO_PASSWORD \\
        --annotation-id ${annotation_id} \\
        --linking-index ${linking_index} \\
        --out-dir ${out_dir} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_dir = "${prefix}"
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/annotation_0.csv
    """
}

process BIOFORMATS2RAWCOMPANION {
    tag "${meta.id}"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/523f11e0352a16c92c05b2baaccfc8e42c7b4d47e4997820478a90f45efa2fbc/data'
        : 'community.wave.seqera.io/library/bioformats2raw:0.9.4--3eec45888b3759e5'}"

    input:
    tuple val(meta), path(master_file), path(master_file_parent)

    output:
    tuple val(meta), path("*.ome.zarr"), emit: ome_zarr
    tuple val("${task.process}"), val('bioformats2rawcompanion'), eval("bioformats2raw --version"), topic: versions, emit: versions_bioformats2rawcompanion

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bioformats2raw \\
        ${args} \\
        ${master_file_parent}/${master_file.name} \\
        --max_workers ${task.cpus} \\
        ${prefix}.ome.zarr \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}.ome.tif
    """
}

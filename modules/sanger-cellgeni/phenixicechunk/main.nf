process PHENIXICECHUNK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c79c0258e345c2e9af8120b1a910a6b0646ca27e9f6f312a085e4fe4f797c809/data'
        : 'community.wave.seqera.io/library/python_procps-ng_pip_phenix-icechunk:af533c50a7c65dc3'}"

    input:
    tuple val(meta), path(image_dir), val(index_name)

    output:
    tuple val(meta), path("*.icechunk"), emit: repo
    tuple val(meta), path("*.summary.json"), emit: summary
    tuple val(meta), path("*.ome.tiff"), emit: ome_tiff, optional: true
    tuple val("${task.process}"), val('phenix-icechunk'), eval("python -c 'import importlib.metadata as m; print(m.version(\"phenix-icechunk\"))'"), emit: versions_phenixicechunk, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The repo holds virtual chunks pointing at the TIFFs under the resolved image_dir,
    // so it stays readable only while the source export stays in place.
    // Per-well OME-TIFFs are emitted when ext.args contains `--ome-tiff .`
    """
    phenix-icechunk \\
        "${image_dir}/${index_name}" \\
        ${prefix}.icechunk \\
        --workers ${task.cpus} \\
        ${args} \\
        > ${prefix}.summary.json
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}.icechunk
    touch ${prefix}.icechunk/repo
    echo '{}' > ${prefix}.summary.json
    touch ${prefix}_A01.ome.tiff
    """
}

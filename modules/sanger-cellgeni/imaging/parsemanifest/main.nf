process IMAGING_PARSEMANIFEST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? '/nfs/cellgeni/singularity/images/stitching_processing.sif'
        : 'biocontainers/YOUR-TOOL-HERE'}"

    input:
    tuple val(meta), path(manifest)

    output:
    tuple val(meta), path("${out_json}"), emit: manifest
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_json = "${prefix}_manifest.json"
    """
    parse_manifest.py --manifest "${manifest}" \\
        --out ${out_json} \\
        ${args} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imaging: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}

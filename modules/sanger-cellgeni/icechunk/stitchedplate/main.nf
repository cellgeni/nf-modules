process ICECHUNK_STITCHEDPLATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/phenix-icechunk:0.3.1'
        : 'quay.io/cellgeni/phenix-icechunk:0.3.1'}"

    input:
    tuple val(meta), path(images, stageAs: 'images/image*/*'), val(layout), path(repos, stageAs: 'sources/source*/*'), val(rounds)

    output:
    tuple val(meta), path("${prefix}.icechunk"), emit: repo
    tuple val(meta), path("${prefix}.summary.json"), emit: summary
    tuple val("${task.process}"), val('phenix-icechunk'), eval("python -c 'import importlib.metadata as m; print(m.version(\"phenix-icechunk\"))'"), emit: versions_phenixicechunk, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_stitched"
    // layout: one [ well, rounds ] map per image, rounds: one round_index per source repo.
    // Only the metadata of the source repos is read, so their source export need not be mounted.
    def image_list = images instanceof List ? images : [images]
    def repo_list = repos instanceof List ? repos : [repos]
    def manifest = groovy.json.JsonOutput.toJson([
        sources: [repo_list, rounds].transpose().collect { repo, round -> [repo: repo.toString(), round_index: round.toString()] },
        images: [image_list, layout].transpose().collect { image, entry -> [path: image.toString(), well: entry.well, rounds: entry.rounds.collect { round -> round.toString() }] },
    ])
    """
    cat <<'END_MANIFEST' > manifest.json
    ${manifest}
    END_MANIFEST

    stitched_plate.py \\
        ${prefix}.icechunk \\
        --manifest manifest.json \\
        --workers ${task.cpus} \\
        ${args} \\
        > ${prefix}.summary.json
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_stitched"
    """
    echo ${args}

    mkdir -p ${prefix}.icechunk
    touch ${prefix}.icechunk/repo
    echo '{}' > ${prefix}.summary.json
    """
}

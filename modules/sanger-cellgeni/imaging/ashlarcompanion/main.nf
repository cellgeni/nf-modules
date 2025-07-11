process IMAGING_ASHLARCOMPANION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ashlar:1.18.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/ashlar:1.18.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(companion_names), path(images)
    path opt_dfp, stageAs: 'dfp*/*'
    path opt_ffp, stageAs: 'ffp*/*'
    val opt_is_plate

    output:
    tuple val(meta), path("${expected_output}"), emit: tif
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_plate = opt_is_plate ? "--plates" : ""
    expected_output = is_plate ? "*/*_${prefix}.ome.tif" : "${prefix}.ome.tif"
    def dfp = opt_dfp ? "--dfp ${opt_dfp}" : ""
    def ffp = opt_ffp ? "--ffp ${opt_ffp}" : ""
    def num_files = companion_names instanceof List ? companion_names.size() : 1
    def opt_dfp_size = opt_dfp instanceof List ? opt_dfp.size() : 1
    def opt_ffp_size = opt_ffp instanceof List ? opt_ffp.size() : 1
    def dfp_validated = opt_dfp_size == 0 || opt_dfp_size == 1 || opt_dfp_size == num_files ? true : false
    def ffp_validated = opt_ffp_size == 0 || opt_ffp_size == 1 || opt_ffp_size == num_files ? true : false

    if (!dfp_validated) {
        error("Please input only zero, one, or N dfp files, where N is the number of input images")
    }
    if (!ffp_validated) {
        error("Please input only zero, one, or N ffp files, where N is the number of input images")
    }

    """

    ashlar \\
        -o ${prefix}.ome.tif \\
        ${companion_names} \\
        ${is_plate} \\
        ${args} \\
        ${dfp} \\
        ${ffp}

    for file in ${expected_output}; do
        sed -i -E 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' "\${file}"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(ashlar --version | sed 's/^.*ashlar //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_plate = opt_is_plate ? "--plates" : ""
    expected_output = is_plate ? "*/*_${prefix}.ome.tif" : "${prefix}.ome.tif"
    def num_files = companion_names instanceof List ? companion_names.size() : 1
    def opt_dfp_size = opt_dfp instanceof List ? opt_dfp.size() : 1
    def opt_ffp_size = opt_ffp instanceof List ? opt_ffp.size() : 1
    def dfp_validated = opt_dfp_size == 0 || opt_dfp_size == 1 || opt_dfp_size == num_files ? true : false
    def ffp_validated = opt_ffp_size == 0 || opt_ffp_size == 1 || opt_ffp_size == num_files ? true : false

    if (!dfp_validated) {
        error("Please input only zero, one, or N dfp files, where N is the number of input images")
    }
    if (!ffp_validated) {
        error("Please input only zero, one, or N ffp files, where N is the number of input images")
    }

    """
    touch ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(ashlar --version | sed 's/^.*ashlar //' )
    END_VERSIONS
    """
}

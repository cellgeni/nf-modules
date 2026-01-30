process SPATIAL_DOT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/dot:latest'
        : 'quay.io/cellgeni/dot:latest'}"

    input:
    tuple val(meta), path(sc_ref), path(sc_spatial)

    output:
    tuple val(meta), path("${prefix}/DOT_weights/DOT_weights_${stem}.csv"), emit: dot_weights
    tuple val(meta), path("${prefix}/DOT_plotdata/DOT_plotting_${stem}.csv"), emit: dot_plotdata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    stem = sc_spatial.baseName
    """
    run_dot.R \\
        ${sc_ref} ${sc_spatial} ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DOT: 0.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    stem = sc_spatial.baseName
    """
    echo ${args}
    
    mkdir -p ${prefix}/DOT_weights
    mkdir -p ${prefix}/DOT_plotdata
    touch ${prefix}/DOT_weights/DOT_weights_${stem}.csv
    touch ${prefix}/DOT_plotdata/DOT_plotting_${stem}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DOT: 0.1.0
    END_VERSIONS
    """
}

process MULTIVIEWSTITCHER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/multiview-stitcher:0.1.62'
        : 'quay.io/cellgeni/multiview-stitcher:0.1.62'}"

    input:
    tuple val(meta), path(repo)

    output:
    tuple val(meta), path("*.ome.zarr"), emit: ome_zarr
    tuple val(meta), path("*_stitched.icechunk"), emit: icechunk
    tuple val(meta), path("*_field_shifts.csv"), emit: shifts
    tuple val(meta), path("*_qc.png"), emit: qc
    tuple val("${task.process}"), val('multiview-stitcher'), eval("python -c 'import importlib.metadata as m; print(m.version(\"multiview-stitcher\"))'"), emit: versions_multiviewstitcher, topic: versions
    tuple val("${task.process}"), val('phenix-icechunk'), eval("python -c 'import importlib.metadata as m; print(m.version(\"phenix-icechunk\"))'"), emit: versions_phenixicechunk, topic: versions

    script:
    def args = task.ext.args ?: ''
    // Outputs are named after the plate in the repo: <plate>_<well>[_z<method>].ome.zarr etc.
    // Stitch a subset of wells with ext.args, e.g. `--wells A4 B3`.
    // itkwasm and matplotlib write caches under HOME, which may not be writable in the container.
    """
    export XDG_CACHE_HOME=\$PWD/.cache
    export MPLCONFIGDIR=\$PWD/.cache/matplotlib

    stitch_wells.py \\
        ${repo} \\
        . \\
        --workers ${task.cpus} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}_A1_zmax.ome.zarr ${prefix}_zmax_stitched.icechunk
    touch ${prefix}_A1_zmax.ome.zarr/zarr.json ${prefix}_zmax_stitched.icechunk/repo
    touch ${prefix}_zmax_field_shifts.csv
    touch ${prefix}_A1_zmax_qc.png
    """
}

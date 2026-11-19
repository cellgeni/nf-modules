process STALIGN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/stalign:2.0.0-spatialdata0.8.0'
        : 'quay.io/cellgeni/stalign:2.0.0-spatialdata0.8.0'}"

    input:
    // source_cells: ordered list of sections (one is fine): cell tables, Xenium output folders or
    // SpatialData .zarr, mixed if needed. Each is staged in its own folder because Xenium outputs
    // all share the name cells.csv.gz
    tuple val(meta), path(source_cells, stageAs: 'sources/???/*'), path(reference_cells, stageAs: 'reference/*'), path(landmarks)

    output:
    tuple val(meta), path("*.zarr"), emit: spatialdata
    tuple val(meta), path("*_aligned_cells.csv.gz"), emit: cells
    tuple val(meta), path("*_transform.npz"), emit: transform
    tuple val(meta), path("*_stalign.json"), emit: summary
    tuple val(meta), path("*_qc.png"), emit: qc, optional: true
    tuple val("${task.process}"), val('stalign'), eval("python -c 'import STalign; print(STalign.__version__)'"), emit: versions_stalign, topic: versions
    tuple val("${task.process}"), val('spatialdata'), eval("python -c 'import importlib.metadata as m; print(m.version(\"spatialdata\"))'"), emit: versions_spatialdata, topic: versions

    script:
    // STalign 2.0.0 is installed from a git commit with its C++ LDDMM engine compiled in the image.
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("STALIGN module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def landmarks_arg = landmarks ? "--landmarks ${landmarks}" : ''
    // Every source section is aligned to the reference and all go into one SpatialData object.
    // Cell tables need x_centroid/y_centroid in micrometers (Xenium cells.csv.gz). SpatialData inputs
    // give their cells from a shapes element (--sdata-shapes, --sdata-cs) and their elements are
    // copied into the output, registered to the reference (--sdata-elements all|none|a,b). Section names
    // (`--names a b c`, default: file names, or section_01.. when these clash), other columns, the
    // method (align_cells or lddmm), LDDMM parameters and the registered affine go in ext.args,
    // e.g. `--names s1 s2 --mixture none`. Sections need not overlap the reference: without landmarks
    // each fit starts by moving the section's mean cell position onto the reference's (--init
    // centroid; --init identity for STalign's default). See stalign_align.py --help.
    """
    export MPLCONFIGDIR=\$PWD/.cache/matplotlib
    export XDG_CACHE_HOME=\$PWD/.cache

    stalign_align.py \\
        --reference ${reference_cells} \\
        --sources ${source_cells} \\
        ${landmarks_arg} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}.zarr
    touch ${prefix}.zarr/zarr.json
    echo "" | gzip > ${prefix}_aligned_cells.csv.gz
    touch ${prefix}_section_01_transform.npz
    touch ${prefix}_stalign.json
    touch ${prefix}_section_01_qc.png
    """
}

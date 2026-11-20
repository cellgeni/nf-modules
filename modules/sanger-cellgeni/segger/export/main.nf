process SEGGER_EXPORT {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/segger:1.0.0'
        : 'quay.io/cellgeni/segger:1.0.0'}"

    input:
    tuple val(meta), path(segmentation), path(source_dir)

    output:
    tuple val(meta), path("${out_dir}"), emit: export_dir
    tuple val(meta), path("${out_dir}/adata.h5ad"), emit: adata, optional: true
    tuple val(meta), path("${out_dir}/cell_boundaries.parquet"), emit: boundaries, optional: true
    tuple val(meta), path("${out_dir}/transcripts.parquet"), emit: transcripts, optional: true
    tuple val(meta), path("${out_dir}/xenium_transcript_assignment.csv"), emit: xenium_transcript_assignment, optional: true
    tuple val(meta), path("${out_dir}/xenium_cell_boundaries.geojson"), emit: xenium_cell_boundaries, optional: true
    tuple val("${task.process}"), val('segger'), eval('python -c "from importlib.metadata import version; print(version(\'segger\'))"'), emit: versions_segger, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger_export"
    """
    segger export \\
        -s ${segmentation} \\
        -i ${source_dir} \\
        -o ${out_dir} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_segger_export"
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/adata.h5ad
    touch ${out_dir}/cell_boundaries.parquet
    touch ${out_dir}/transcripts.parquet
    touch ${out_dir}/xenium_transcript_assignment.csv
    printf '%s\\n' '{"type":"FeatureCollection","features":[]}' > ${out_dir}/xenium_cell_boundaries.geojson
    """
}

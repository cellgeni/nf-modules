process BAYSOR_RUN {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/baysor:cpp-0.8.2"

    input:
    tuple val(meta), path(coordinates), path(prior_segmentation)
    path config

    output:
    tuple val(meta), path("${out_dir}"                                                                      , type: 'dir'), emit: outdir
    tuple val(meta), path("${out_dir}/{segmentation.csv,molecules.parquet}")                                               , emit: molecules
    tuple val(meta), path("${out_dir}/{segmentation_cell_stats.csv,cells.parquet}")                                        , emit: cell_stats
    tuple val(meta), path("${out_dir}/{segmentation_counts.loom,segmentation_counts.tsv,feature_matrix.h5}")               , emit: counts
    tuple val(meta), path("${out_dir}/{segmentation_polygons_2d.json,cell_boundaries.parquet}")             , optional: true, emit: polygons_2d
    tuple val(meta), path("${out_dir}/{segmentation_polygons_3d.json,cell_boundaries_3d.parquet}")          , optional: true, emit: polygons_3d
    tuple val(meta), path("${out_dir}/segmentation_plot.html")                                              , optional: true, emit: plot
    tuple val(meta), path("${out_dir}/diagnostic_report.html")                                              , optional: true, emit: report
    tuple val(meta), path("${out_dir}/{segmentation_params.dump.toml,run_params.toml}")                                    , emit: params
    tuple val(meta), path("${out_dir}/{segmentation_log.log,run.log}")                                                     , emit: log
    tuple val("${task.process}"), val('baysor'), val('0.8.2'), topic: versions, emit: versions_baysor_run

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Baysor takes the prior segmentation as a second positional argument. It is either a staged
    // file (image mask / boundary table) or the literal ':column_name' form, which cannot be a
    // Nextflow path input and is therefore passed through task.ext.prior instead.
    def prior = prior_segmentation ? "${prior_segmentation}" : (task.ext.prior ?: '')
    def config_arg = config ? "--config ${config}" : ''
    out_dir = "${prefix}_baysor"
    """
    export OMP_NUM_THREADS=${task.cpus}

    baysor run \\
        ${config_arg} \\
        --output ${out_dir} \\
        ${args} \\
        ${coordinates} \\
        ${prior}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_baysor"
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/segmentation.csv
    touch ${out_dir}/segmentation_cell_stats.csv
    touch ${out_dir}/segmentation_counts.loom
    touch ${out_dir}/segmentation_polygons_2d.json
    touch ${out_dir}/segmentation_params.dump.toml
    touch ${out_dir}/segmentation_log.log
    """
}

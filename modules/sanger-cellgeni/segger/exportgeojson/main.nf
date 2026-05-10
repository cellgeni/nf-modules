process SEGGER_EXPORTGEOJSON {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/cellgeni/segger:v2"

    input:
    tuple val(meta), path(zarr_zip)

    output:
    tuple val(meta), path("${prefix}.geojson"), emit: geojson
    tuple val("${task.process}"), val('export_zarr_polygons_to_geojson'), eval('export_zarr_polygons_to_geojson.py --version'), emit: versions_export_zarr_polygons_to_geojson, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export_zarr_polygons_to_geojson.py \\
        "${zarr_zip}" \\
        "${prefix}.geojson" \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\n' '{"type":"FeatureCollection","features":[]}' > "${prefix}.geojson"
    """
}

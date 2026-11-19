process GEOJSON_FROMXENIUMZARR {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'quay.io/cellgeni/spatialdata:0.9.0-spatialdata-io0.7.2'
        : 'quay.io/cellgeni/spatialdata:0.9.0-spatialdata-io0.7.2'}"

    input:
    tuple val(meta), path(zarr_zip)

    output:
    tuple val(meta), path("${prefix}.geojson"), emit: geojson
    tuple val("${task.process}"), val('export_zarr_polygons_to_geojson'), eval('export_zarr_polygons_to_geojson.py --version'), emit: versions_export_zarr_polygons_to_geojson, topic: versions

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

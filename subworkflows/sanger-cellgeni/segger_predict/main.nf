include { SEGGER_SEGMENT       } from '../../../modules/sanger-cellgeni/segger/segment/main'
include { SEGGER_EXPORT        } from '../../../modules/sanger-cellgeni/segger/export/main'
include { SEGGER_EXPORTGEOJSON } from '../../../modules/sanger-cellgeni/segger/exportgeojson/main'

workflow SEGGER_PREDICT {
    take:
    ch_input // channel: [ val(meta), path(input_dir) ]

    main:
    def ch_versions = channel.empty()

    SEGGER_SEGMENT(ch_input)
    ch_versions = ch_versions.mix(SEGGER_SEGMENT.out.versions.first())

    // Join segmentation result with the original input_dir to use as source_dir for export
    def ch_export = SEGGER_SEGMENT.out.segmentation.join(ch_input, by: 0)

    SEGGER_EXPORT(ch_export)
    ch_versions = ch_versions.mix(SEGGER_EXPORT.out.versions.first())

    // SEGGER_EXPORT writes Xenium-compatible output into export_dir; export the cell polygons to GeoJSON.
    def ch_polygons_zarr = SEGGER_EXPORT.out.export_dir.map { meta, export_dir -> [meta, file("${export_dir}/seg_cells.zarr.zip")] }

    SEGGER_EXPORTGEOJSON(ch_polygons_zarr)
    ch_versions = ch_versions.mix(SEGGER_EXPORTGEOJSON.out.versions_export_zarr_polygons_to_geojson.first())
    def ch_segmentation = SEGGER_SEGMENT.out.segmentation
    def ch_export_dir = SEGGER_EXPORT.out.export_dir
    def ch_geojson = SEGGER_EXPORTGEOJSON.out.geojson

    emit:
    segmentation = ch_segmentation // channel: [ val(meta), path(parquet) ]
    export_dir   = ch_export_dir // channel: [ val(meta), path(export_dir) ]
    geojson      = ch_geojson // channel: [ val(meta), path(geojson) ]
    versions     = ch_versions // channel: [ versions.yml ]
}

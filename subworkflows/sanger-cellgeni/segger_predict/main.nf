include { SEGGER_SEGMENT } from '../../../modules/sanger-cellgeni/segger/segment/main'
include { SEGGER_EXPORT  } from '../../../modules/sanger-cellgeni/segger/export/main'

workflow SEGGER_PREDICT {
    take:
    ch_input // channel: [ val(meta), path(input_dir) ]

    main:
    def ch_versions = channel.empty()

    SEGGER_SEGMENT(ch_input)
    ch_versions = ch_versions.mix(SEGGER_SEGMENT.out.versions_segger.first())

    // Join segmentation result with the original input_dir to use as source_dir for export
    def ch_export = SEGGER_SEGMENT.out.segmentation.join(ch_input, by: 0)

    // The GeoJSON and transcript assignment are only written when SEGGER_EXPORT's
    // ext.args include the `xeniumranger` element (Xenium input only).
    SEGGER_EXPORT(ch_export)
    ch_versions = ch_versions.mix(SEGGER_EXPORT.out.versions_segger.first())

    emit:
    segmentation          = SEGGER_SEGMENT.out.segmentation                // channel: [ val(meta), path(parquet) ]
    export_dir            = SEGGER_EXPORT.out.export_dir                   // channel: [ val(meta), path(export_dir) ]
    geojson               = SEGGER_EXPORT.out.xenium_cell_boundaries       // channel: [ val(meta), path(geojson) ]
    transcript_assignment = SEGGER_EXPORT.out.xenium_transcript_assignment // channel: [ val(meta), path(csv) ]
    versions              = ch_versions                                     // channel: [ versions ]
}
